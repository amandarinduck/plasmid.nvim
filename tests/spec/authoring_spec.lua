local H = require("helpers")
local pio = require("plasmid.io")
local tx = require("plasmid.model.transaction")
local author = require("plasmid.author")
local sel = require("plasmid.select")
local sg = require("plasmid.io.snapgene")

require("plasmid").setup({})
local buffer = require("plasmid.buffer")

local T = { _order = {} }
local function case(name, fn) T[name] = fn; table.insert(T._order, name) end

local tmp = vim.fn.tempname()
vim.fn.mkdir(tmp, "p")
local function copy(rel, name)
  local dst = tmp .. "/" .. (name or vim.fn.fnamemodify(rel, ":t"))
  local fh = assert(io.open(dst, "wb")); fh:write(H.read(H.fixture(rel))); fh:close()
  return dst
end
local function open(path)
  vim.cmd("silent edit " .. vim.fn.fnameescape(path))
  local buf = vim.api.nvim_get_current_buf()
  return buf, buffer.get(buf), buffer.record(buffer.get(buf))
end
local function keys(s) vim.api.nvim_feedkeys(vim.api.nvim_replace_termcodes(s, true, false, true), "x", false) end

-- model snapshot for exact-state comparisons
local function snapshot(rec)
  local out = { rec.topology }
  for _, f in ipairs(rec.features) do
    local q = {}
    for _, x in ipairs(f.qualifiers) do q[#q + 1] = x.key .. "=" .. tostring(x.value) end
    out[#out + 1] = table.concat({ f.id, f.type, f.location_text, table.concat(q, ";") }, "|")
  end
  return out
end

vim.o.columns, vim.o.lines = 120, 50
vim.o.hidden = true

case("location_for_selections", function()
  local rec = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  H.eq("100..200", author.location_for(rec, sel.make(rec, 100, 200), 1))
  H.eq("complement(100..200)", author.location_for(rec, sel.make(rec, 100, 200), -1))
  H.eq("join(4300..4361,1..50)", author.location_for(rec, sel.make(rec, 4300, 50), 1))
  H.eq("complement(join(4300..4361,1..50))", author.location_for(rec, sel.make(rec, 4300, 50), -1))
  H.eq("77", author.location_for(rec, sel.make(rec, 77, 77), 1))
end)

case("editor_text_roundtrips_every_corpus_feature", function()
  local n = 0
  for _, path in ipairs({ "ncbi/J01749.gb", "ncbi/U01317.gb", "ncbi/NC_001422.gb", "synthetic/grammar.gb",
    "snapgene/pFA-KanMX4.dna", "snapgene/linebreak_in_qualifier_text.dna" }) do
    local rec = assert(pio.read(H.fixture(path))).records[1]
    for _, f in ipairs(rec.features) do
      local lines = author.editor_text(rec, f, "x")
      local g, err = author.parse_editor(rec, lines)
      H.ok(g, path .. " " .. f.location_text .. ": " .. tostring(err))
      H.eq(f.type, g.type)
      H.eq(f.location_text, g.location_text)
      H.eq(#f.qualifiers, #g.qualifiers, path .. " " .. f.location_text)
      for i, q in ipairs(f.qualifiers) do
        H.eq({ q.key, q.value }, { g.qualifiers[i].key, g.qualifiers[i].value }, path .. " " .. f.location_text)
      end
      n = n + 1
    end
  end
  H.ok(n > 100, "features checked: " .. n)
end)

case("editor_parsing_is_forgiving_and_validates", function()
  local rec = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local f = assert(author.parse_editor(rec, {
    "# comment", "", "  promoter   complement(join(10..20,", "30..40))", '/label="my ""best"" promoter"',
    "   /note=one line", "/pseudo",
  }))
  H.eq("promoter", f.type)
  H.eq("complement(join(10..20,30..40))", f.location_text)
  H.eq('my "best" promoter', f:get("label"))
  H.eq(true, f:get("pseudo"))
  local _, e1 = author.parse_editor(rec, { "CDS join(1..10" })
  H.ok(e1 and e1:match("invalid location"), tostring(e1))
  local _, e2 = author.parse_editor(rec, { "CDS 4000..5000" })
  H.ok(e2 and e2:match("outside 1..4361"), tostring(e2))
  local _, e3 = author.parse_editor(rec, { "# only comments" })
  H.ok(e3 and e3:match("empty"), tostring(e3))
end)

case("undo_redo_restore_exact_state_property", function()
  math.randomseed(11)
  local rec = assert(pio.read(H.fixture("ncbi/J01749.gb"))).records[1]
  local base = snapshot(rec)
  local states = { base }
  for step = 1, 60 do
    local r = math.random()
    if r < 0.4 or #rec.features < 3 then
      local a = math.random(1, 4000)
      local f = author.make_feature(rec, { type = "misc_feature", label = "f" .. step,
        sel = sel.make(rec, a, a + math.random(0, 300)), strand = math.random() < 0.5 and 1 or -1 })
      tx.commit(rec, { op = "add", feature = f })
    elseif r < 0.7 then
      local f = rec.features[math.random(1, #rec.features)]
      tx.commit(rec, { op = "remove", id = f.id })
    elseif r < 0.9 then
      local f = rec.features[math.random(1, #rec.features)]
      local g = assert(author.parse_editor(rec, { "misc_feature 5..25", '/label="edited ' .. step .. '"' }))
      tx.commit(rec, { op = "replace", id = f.id, feature = g })
    else
      tx.commit(rec, { op = "topology", value = rec.topology == "linear" and "circular" or "linear" })
    end
    states[#states + 1] = snapshot(rec)
  end
  for i = #states - 1, 1, -1 do
    H.ok(tx.undo(rec))
    H.eq(states[i], snapshot(rec), "after undo to state " .. i)
  end
  H.eq(nil, tx.undo(rec))
  for i = 2, #states do
    H.ok(tx.redo(rec))
    H.eq(states[i], snapshot(rec), "after redo to state " .. i)
  end
end)

case("genbank_save_keeps_untouched_features_verbatim", function()
  local path = copy("ncbi/J01749.gb", "edit.gb")
  local original = H.read(path)
  local res = assert(pio.read(path))
  local rec = res.records[1]
  local f = author.make_feature(rec, { type = "misc_feature", label = "my site", sel = sel.make(rec, 1500, 1520), strand = -1 })
  tx.commit(rec, { op = "add", feature = f })
  local out = assert(pio.format(res.records, "genbank", res.file))
  -- removing exactly the new feature's canonical lines restores the original bytes
  local added = table.concat(require("plasmid.io.genbank").format_feature(f), "\n") .. "\n"
  local s, e = out:find(added, 1, true)
  H.ok(s, "new feature lines not found")
  H.eq(original, out:sub(1, s - 1) .. out:sub(e + 1))
  local back = assert(pio.parse(out, path)).records[1]
  local g
  for _, x in ipairs(back.features) do if x:label() == "my site" then g = x end end
  H.eq("complement(1500..1520)", g.location_text)
  -- undo, and the file is byte-identical again
  tx.undo(rec)
  H.eq(original, assert(pio.format(res.records, "genbank", res.file)))
end)

local function packets_except(bytes, skip)
  local out = {}
  for _, p in ipairs(assert(sg.packets(bytes))) do
    if p.type ~= skip then out[#out + 1] = p.type .. ":" .. p.data end
  end
  return out
end

case("snapgene_save_add_edit_delete", function()
  for _, rel in ipairs({ "snapgene/pFA-KanMX4.dna", "snapgene/sample-d.dna", "snapgene/sgffp_origin_spanning_features.dna" }) do
    local path = copy(rel, "w_" .. vim.fn.fnamemodify(rel, ":t"))
    local original = H.read(path)
    local res = assert(pio.read(path))
    local rec = res.records[1]
    local L = rec:len()
    local nprimers = 0
    for _, f in ipairs(rec.features) do if f.origin == "primer" then nprimers = nprimers + 1 end end
    -- add a reverse feature, an origin-wrapping one, edit one, delete one
    tx.commit(rec, { op = "add", feature = author.make_feature(rec, { type = "protein_bind", label = "LacO <test> & \"q\"",
      sel = sel.make(rec, 20, 40), strand = -1, color = "#33aa77" }) })
    if rec:is_circular() then
      tx.commit(rec, { op = "add", feature = author.make_feature(rec, { type = "misc_feature", label = "wrap",
        sel = sel.make(rec, L - 9, 12), strand = 1 }) })
    end
    local victim, edited
    for _, f in ipairs(rec.features) do
      if f.origin ~= "primer" and f.raw_xml then
        if not victim then victim = f elseif not edited then edited = f end
      end
    end
    tx.commit(rec, { op = "remove", id = victim.id })
    if edited then
      local g = assert(author.parse_editor(rec, { "CDS 5..37", '/label="renamed"', "/codon_start=1" }))
      tx.commit(rec, { op = "replace", id = edited.id, feature = g })
    end
    local out = assert(pio.format(res.records, "snapgene", res.file))
    -- every packet other than Features is byte-identical
    H.eq(packets_except(original, 0x0A), packets_except(out, 0x0A), rel)
    -- re-reading gives back the model (primers come from their own packet)
    local back = assert(pio.parse(out, path)).records[1]
    local function feats(r)
      local t = {}
      for _, f in ipairs(r.features) do
        if f.origin ~= "primer" then
          local q = {}
          for _, x in ipairs(f.qualifiers) do
            if x.key ~= "ApEinfo_fwdcolor" and x.key ~= "ApEinfo_revcolor" then q[#q + 1] = x.key .. "=" .. tostring(x.value) end
          end
          table.sort(q)
          t[#t + 1] = f.type .. "|" .. f.location_text .. "|" .. f:label() .. "|" .. table.concat(q, ";")
        end
      end
      return t
    end
    H.eq(feats(rec), feats(back), rel)
    local count = 0
    for _, f in ipairs(back.features) do if f.origin == "primer" then count = count + 1 end end
    H.eq(nprimers, count, rel .. " primers")
    for _, f in ipairs(back.features) do
      if f:label() == "LacO <test> & \"q\"" then
        H.eq(-1, f.strand)
        H.eq("#33aa77", f.color)
      end
    end
    -- persist for the Biopython cross-check fixture set
    local keep = H.root .. "/tests/fixtures/written/" .. vim.fn.fnamemodify(rel, ":t")
    vim.fn.mkdir(vim.fn.fnamemodify(keep, ":h"), "p")
    if vim.env.PLASMID_REGEN_WRITTEN then
      local fh = assert(io.open(keep, "wb")); fh:write(out); fh:close()
    elseif vim.fn.filereadable(keep) == 1 then
      H.eq(H.read(keep), out, rel .. ": writer output changed; regenerate fixtures and re-run the Biopython check")
    end
    -- undo everything: Features packet regenerated from verbatim XML equals the original
    while tx.undo(rec) do end
    H.eq(original, assert(pio.format(res.records, "snapgene", res.file)), rel .. " after undo")
  end
end)

case("snapgene_refuses_edited_primers", function()
  local res = assert(pio.read(H.fixture("snapgene/sample-hybridization-params.dna")))
  local rec = res.records[1]
  local primer
  for _, f in ipairs(rec.features) do if f.origin == "primer" then primer = f break end end
  H.ok(primer, "fixture has primers")
  tx.commit(rec, { op = "remove", id = primer.id })
  local out, err = pio.format(res.records, "snapgene", res.file)
  H.eq(nil, out)
  H.ok(err:match("primer"), err)
  H.ok(pio.format(res.records, "genbank", res.file, "snapgene"))
end)

case("quick_add_undo_redo_in_buffer", function()
  local buf, st, rec = open(copy("ncbi/J01749.gb", "ui_add.gb"))
  local n0 = #rec.features
  vim.cmd("DnaSelect 2000 2100")
  local answers = { "my feature", "promoter", "-" }
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb(table.remove(answers, 1)) end
  keys("ga")
  vim.ui.input = orig
  H.eq(n0 + 1, #rec.features)
  local f
  for _, x in ipairs(rec.features) do if x:label() == "my feature" then f = x end end
  H.eq("promoter", f.type)
  H.eq("complement(2000..2100)", f.location_text)
  H.ok(f.pinned)
  H.eq(true, vim.bo[buf].modified)
  H.ok(buffer.winbar():find("unsaved", 1, true))
  -- the pinned feature is drawn even when tracks are limited to 1
  vim.cmd("DnaTracks 1")
  local drawn = false
  for _, info in pairs(st.layout.rows) do
    for _, it in ipairs(info.items or {}) do if it.f == f then drawn = true end end
  end
  H.ok(drawn, "pinned feature hidden by track limit")
  vim.cmd("DnaTracks")
  keys("u")
  H.eq(n0, #rec.features)
  H.eq(false, vim.bo[buf].modified)
  keys("<C-r>")
  H.eq(n0 + 1, #rec.features)
  vim.cmd("silent write")
  H.eq(false, vim.bo[buf].modified)
  local back = assert(pio.read(st.path)).records[1]
  H.eq(n0 + 1, #back.features)
end)

case("visual_ga_and_editor_flow", function()
  local buf, st, rec = open(copy("ncbi/J01749.gb", "ui_edit.gb"))
  local win = vim.api.nvim_get_current_win()
  -- visual range -> ga
  buffer.goto_bp(win, 300)
  local answers = { "visual feat", "misc_feature", "+" }
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb(table.remove(answers, 1)) end
  keys("v<Cmd>lua require('plasmid.buffer').goto_bp(0, 340)<CR>ga")
  vim.ui.input = orig
  local f
  for _, x in ipairs(rec.features) do if x:label() == "visual feat" then f = x end end
  H.eq("300..340", f.location_text)
  -- edit it through the editor buffer
  require("plasmid.ui").goto_feature(buf, f.id)
  keys("gE")
  local ebuf = vim.api.nvim_get_current_buf()
  H.ok(vim.api.nvim_buf_get_name(ebuf):match("^plasmid://feature/"))
  vim.api.nvim_buf_set_lines(ebuf, 5, -1, false, { "terminator complement(300..345)", '/label="edited feat"', "/note=hello" })
  vim.cmd("silent write")
  H.eq(buf, vim.api.nvim_get_current_buf(), "editor window should close on success")
  local g = rec:feature_by_id(f.id)
  H.eq({ "terminator", "complement(300..345)", "edited feat", "hello" }, { g.type, g.location_text, g:label(), g:get("note") })
  -- a bad edit keeps the editor open and changes nothing
  require("plasmid.ui").goto_feature(buf, f.id)
  keys("gE")
  local ebuf2 = vim.api.nvim_get_current_buf()
  vim.api.nvim_buf_set_lines(ebuf2, 5, -1, false, { "terminator 1..99999" })
  local before = #H.messages
  pcall(vim.cmd, "silent write")
  local reported = false
  for i = before + 1, #H.messages do
    local m = H.messages[i]
    if m.level == vim.log.levels.ERROR and tostring(m.msg):match("outside 1..4361") then reported = true end
  end
  H.ok(reported, "invalid location was not reported as an error")
  H.eq(ebuf2, vim.api.nvim_get_current_buf())
  H.eq("complement(300..345)", rec:feature_by_id(f.id).location_text)
  vim.cmd("bwipeout!")
  -- delete with confirmation stubbed, then undo
  require("plasmid.ui").goto_feature(buf, f.id)
  local osel = vim.ui.select
  vim.ui.select = function(items, _, cb) cb(items[1], 1) end
  vim.cmd("DnaFeatureDelete")
  vim.ui.select = osel
  H.eq(nil, rec:feature_by_id(f.id))
  vim.cmd("DnaUndo")
  H.ok(rec:feature_by_id(f.id))
end)

case("personal_library", function()
  require("plasmid.config").options.library_path = tmp .. "/lib/library.gb"
  local buf, st, rec = open(copy("ncbi/J01749.gb", "lib.gb"))
  local bla
  for _, f in ipairs(rec.features) do if f.type == "CDS" and f:get("gene") == "bla" then bla = f end end
  require("plasmid.ui").goto_feature(buf, bla.id)
  vim.cmd("DnaFeatureSave")
  vim.cmd("DnaFeatureSave") -- duplicate is skipped
  local lib = assert(pio.read(tmp .. "/lib/library.gb"))
  H.eq(1, #lib.records)
  H.eq(rec:extract(bla), lib.records[1].seq:tostring())
  H.eq("CDS", lib.records[1].features[1].type)
  H.eq("1..861", lib.records[1].features[1].location_text)
  -- save a plain selection
  vim.cmd("DnaSelectClear")
  vim.cmd("DnaSelect 10 30")
  buffer.goto_bp(0, 2000) -- cursor away from features
  local answers = { "spacer", "misc_feature" }
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb(table.remove(answers, 1)) end
  local orig_fuc = sel.feature_under_cursor
  sel.feature_under_cursor = function() return nil end
  local orig_ci = buffer.cursor_info
  buffer.cursor_info = function() return 2000, nil end
  vim.cmd("DnaFeatureSave")
  vim.ui.input, sel.feature_under_cursor, buffer.cursor_info = orig, orig_fuc, orig_ci
  lib = assert(pio.read(tmp .. "/lib/library.gb"))
  H.eq(2, #lib.records)
  H.eq(rec.seq:sub(10, 30), lib.records[2].seq:tostring())
  require("plasmid.config").options.library_path = nil
end)

case("fasta_save_warns_about_features", function()
  local buf, st, rec = open(copy("ncbi/J01749.gb", "warn.gb"))
  vim.cmd("DnaSelect 100 200")
  local answers = { "x", "misc_feature", "+" }
  local orig = vim.ui.input
  vim.ui.input = function(_, cb) cb(table.remove(answers, 1)) end
  keys("ga")
  vim.ui.input = orig
  local msgs = {}
  local on = vim.notify
  vim.notify = function(m) msgs[#msgs + 1] = m end
  vim.cmd("silent write " .. vim.fn.fnameescape(tmp .. "/warn.fasta"))
  vim.notify = on
  H.ok(table.concat(msgs, "\n"):match("FASTA cannot store features"), table.concat(msgs, "\n"))
end)

return T
