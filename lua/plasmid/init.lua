-- plasmid.nvim: DNA sequence viewing and editing for Neovim.
local config = require("plasmid.config")

local M = {}

local initialized = false
local group

local function cmd(name, fn, opts)
  vim.api.nvim_create_user_command(name, fn, opts or {})
end

local function cur()
  local buffer = require("plasmid.buffer")
  local buf = vim.api.nvim_get_current_buf()
  local st = buffer.get(buf)
  if not st then
    require("plasmid.util").notify("plasmid: current buffer is not a sequence file", vim.log.levels.WARN)
    return nil
  end
  return buf, st, buffer.record(st)
end

local function register_commands()
  local buffer = require("plasmid.buffer")
  local ui = require("plasmid.ui")

  cmd("DnaOpen", function(a)
    vim.cmd.edit(vim.fn.fnameescape(a.args))
  end, { nargs = 1, complete = "file", desc = "Open a FASTA/GenBank file" })

  cmd("DnaView", function(a)
    local buf = cur()
    if not buf then return end
    local v = a.args ~= "" and a.args or "linear"
    if v ~= "linear" and v ~= "circular" then
      return require("plasmid.util").notify("plasmid: view must be linear or circular", vim.log.levels.ERROR)
    end
    buffer.set_view(buf, v, a.bang)
  end, {
    nargs = "?", bang = true, desc = "Switch view (! forces circular on linear records)",
    complete = function() return { "linear", "circular" } end,
  })

  cmd("DnaToggle", function()
    local buf = cur()
    if buf then buffer.toggle_view(buf) end
  end, { desc = "Toggle linear/circular view" })

  cmd("DnaGoto", function(a)
    local buf = cur()
    if not buf then return end
    local bp = tonumber(a.args)
    if not bp then return require("plasmid.util").notify("plasmid: :DnaGoto {bp}", vim.log.levels.ERROR) end
    buffer.goto_bp(vim.api.nvim_get_current_win(), bp)
  end, { nargs = 1, desc = "Go to base pair" })

  cmd("DnaSelect", function(a)
    local buf, st, rec = cur()
    if not buf then return end
    local s = require("plasmid.select")
    local x, y = a.args:match("^%s*(%d+)%s*[%.%s,%-]+%s*(%d+)%s*$")
    if x then return s.set(buf, s.make(rec, tonumber(x), tonumber(y))) end
    local want = a.args:lower()
    for _, f in ipairs(rec.features) do
      if f:label(config.options.label_qualifiers):lower() == want then
        require("plasmid.ui").goto_feature(buf, f.id)
        return s.set(buf, s.for_feature(rec, f))
      end
    end
    require("plasmid.util").notify("plasmid: :DnaSelect {start} {stop} | {feature label}", vim.log.levels.ERROR)
  end, { nargs = "+", desc = "Select a range or a feature by label",
    complete = function()
      local st = require("plasmid.buffer").get(0)
      if not st then return {} end
      local out = {}
      for _, f in ipairs(require("plasmid.buffer").record(st).features) do out[#out + 1] = f:label() end
      return out
    end })
  cmd("DnaYank", function(a)
    local buf = cur()
    if not buf then return end
    local s = require("plasmid.select")
    local key, reg = nil, a.fargs[2]
    if a.fargs[1] then
      key = s.strand_key(a.fargs[1])
      if not key or (reg and #reg ~= 1) then
        return require("plasmid.util").notify("plasmid: :DnaYank [top|top35|bottom|bottom35] [register]", vim.log.levels.ERROR)
      end
    end
    local sel = s.get(buf)
    if not sel then
      return require("plasmid.util").notify("plasmid: nothing selected (vaf, visual mode, or :DnaSelect)", vim.log.levels.WARN)
    end
    s.yank_strand(buf, sel, key, reg)
  end, { nargs = "*", desc = "Yank the selection: top/bottom strand, either direction (no argument: choose)",
    complete = function() return { "top", "top35", "bottom", "bottom35" } end })
  cmd("DnaSelectClear", function()
    local buf = cur()
    if buf then require("plasmid.select").set(buf, nil) end
  end, { desc = "Clear the selection" })

  cmd("DnaOrfs", function(a)
    local buf = cur()
    if not buf then return end
    local orfs = require("plasmid.orfs")
    if a.args == "off" then
      local st = require("plasmid.buffer").get(buf)
      if st.orfs then orfs.toggle(buf) end
      return
    end
    orfs.toggle(buf, tonumber(a.args))
  end, { nargs = "?", desc = "Toggle ORFs; :DnaOrfs {min_aa} shows ORFs of at least that length; :DnaOrfs off",
    complete = function() return { "off", "50", "75", "100", "150" } end })
  cmd("DnaOrfAdd", function()
    if cur() then require("plasmid.orfs").promote_under_cursor() end
  end, { desc = "Turn the ORF under the cursor into a CDS feature" })

  cmd("DnaFeatureDbImport", function(a)
    local args = a.fargs
    local n, name, note = require("plasmid.featuredb").import(args[1], args[2], args[3])
    if not n then return require("plasmid.util").notify("plasmid: import failed: " .. tostring(name), vim.log.levels.ERROR) end
    require("plasmid.util").notify(("plasmid: imported %d reference features as '%s'%s"):format(n, name,
      note and (" (" .. note .. ")") or ""))
  end, { nargs = "+", complete = "file",
    desc = "Import a feature database: {pLannotate BLAST_dbs dir | FASTA [name] [descriptions.csv] | GenBank/.dna file or folder}" })
  cmd("DnaFeatureDbList", function()
    local fdb = require("plasmid.featuredb")
    local lines = { "Feature databases in " .. vim.fn.fnamemodify(fdb.dir(), ":~"), "" }
    for _, n in ipairs(fdb.list()) do lines[#lines + 1] = ("  %-20s %5d features"):format(n, #(fdb.load(n) or {})) end
    if #lines == 2 then lines[#lines + 1] = "  (none; :DnaFeatureDbImport {path})" end
    lines[#lines + 1] = ""
    lines[#lines + 1] = "Personal library: " .. vim.fn.fnamemodify(require("plasmid.author").library_path(), ":~")
    vim.lsp.util.open_floating_preview(lines, "", { border = "rounded", focus_id = "plasmid_db" })
  end, { desc = "List feature databases" })
  cmd("DnaDetect", function(a)
    local buf = cur()
    if not buf then return end
    if a.fargs[1] == "off" then return require("plasmid.detectui").clear(buf) end
    require("plasmid.detectui").run(buf, a.fargs)
  end, { nargs = "*", desc = "Detect known features (optionally limited to named databases); :DnaDetect off",
    complete = function()
      local out = { "off" }
      vim.list_extend(out, require("plasmid.featuredb").list())
      return out
    end })
  cmd("DnaDetectAccept", function(a)
    local buf = cur()
    if not buf then return end
    local d = require("plasmid.detectui")
    if a.bang then return d.accept_all(buf, a.args == "fragments") end
    local _, fid = require("plasmid.buffer").cursor_info(vim.api.nvim_get_current_win())
    if not (fid and d.accept(buf, fid)) then
      require("plasmid.util").notify("plasmid: put the cursor on a proposed feature (or :DnaDetectAccept! for all)",
        vim.log.levels.WARN)
    end
  end, { bang = true, nargs = "?", desc = "Accept the proposal under the cursor (! = all full-length; ! fragments = all)" })
  cmd("DnaDetectReview", function() if cur() then require("plasmid.detectui").review(vim.api.nvim_get_current_buf()) end end,
    { desc = "Pick proposed features to accept" })

  local author = require("plasmid.author")
  cmd("DnaFeatureAdd", function(a)
    local buf, st, rec = cur()
    if not buf then return end
    local s = require("plasmid.select")
    local x, y = a.args:match("^%s*(%d+)%s*[%.%s,%-]+%s*(%d+)%s*$")
    if x then s.set(buf, s.make(rec, tonumber(x), tonumber(y))) end
    if a.bang or not s.get(buf) then return author.add_with_editor() end
    author.quick_add(false)
  end, { nargs = "?", bang = true,
    desc = "Add a feature over the selection or {start}..{stop} (! opens the full editor)" })
  cmd("DnaFeatureEdit", function() if cur() then author.edit_under_cursor() end end,
    { desc = "Edit the feature under the cursor" })
  cmd("DnaFeatureDelete", function(a) if cur() then author.delete(a.bang) end end,
    { bang = true, desc = "Delete the feature under the cursor (! skips confirmation)" })
  cmd("DnaFeatureSave", function() if cur() then author.save_to_library() end end,
    { desc = "Save the feature under the cursor (or selection) to your library" })
  cmd("DnaLibrary", function() vim.cmd.edit(vim.fn.fnameescape(author.library_path())) end,
    { desc = "Open your personal feature library" })
  cmd("DnaUndo", function() if cur() then author.undo() end end, { desc = "Undo last feature change" })
  cmd("DnaRedo", function() if cur() then author.redo() end end, { desc = "Redo" })

  -- Sequence editing ------------------------------------------------------
  -- A range argument is "s..e" (also "s-e", "s e"); without one, the
  -- selection. Returns s, e, rest-of-args or nil.
  local function range_args(buf, args, need_rest)
    local s, e, rest = args:match("^%s*(%d+)%s*%.%.%s*(%d+)%s*(.-)%s*$")
    if not s then s, e, rest = args:match("^%s*(%d+)%s*%-%s*(%d+)%s*(.-)%s*$") end
    if not s and not need_rest then s, e, rest = args:match("^%s*(%d+)%s+(%d+)%s*(.-)%s*$") end
    if s then return tonumber(s), tonumber(e), rest end
    local sel = require("plasmid.select").get(buf)
    if sel then return sel.start, sel.stop, (args:gsub("^%s+", ""):gsub("%s+$", "")) end
    return nil
  end
  local function in_range(rec, ...)
    for _, x in ipairs({ ... }) do
      if x < 1 or x > rec:len() then
        require("plasmid.util").notify(("plasmid: position %d outside 1..%d"):format(x, rec:len()), vim.log.levels.ERROR)
        return false
      end
    end
    return true
  end
  local se = function() return require("plasmid.seqedit") end
  cmd("DnaInsert", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    local pos, text = a.args:match("^%s*(%d+)%s+(.+)$")
    pos = tonumber(pos)
    if not pos or pos < 1 or pos > rec:len() + 1 then
      return require("plasmid.util").notify(("plasmid: :DnaInsert {pos} {residues}  (inserted before pos; 1..%d)"):format(rec:len() + 1),
        vim.log.levels.ERROR)
    end
    se().insert(buf, pos, text, { force = a.bang })
  end, { nargs = "+", bang = true, desc = "Insert residues before a position" })
  cmd("DnaDelete", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    local s, e = range_args(buf, a.args)
    if not s then return require("plasmid.util").notify("plasmid: :DnaDelete[!] [{start}..{stop}] (default: the selection)", vim.log.levels.ERROR) end
    if not in_range(rec, s, e) then return end
    se().delete(buf, s, e, { force = a.bang })
  end, { nargs = "*", bang = true, desc = "Delete a range or the selection (! overrides cut_policy=refuse)" })
  cmd("DnaReplace", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    local s, e, text = range_args(buf, a.args, true)
    if not s or not text or text == "" then
      return require("plasmid.util").notify("plasmid: :DnaReplace[!] [{start}..{stop}] {residues} (default range: the selection)", vim.log.levels.ERROR)
    end
    if not in_range(rec, s, e) then return end
    se().replace(buf, s, e, text, { force = a.bang })
  end, { nargs = "+", bang = true, desc = "Replace a range or the selection by residues" })
  cmd("DnaRevComp", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    if a.args:match("^%s*all%s*$") then return se().revcomp(buf, nil, nil, { force = a.bang }) end
    local s, e = range_args(buf, a.args)
    if not s then return require("plasmid.util").notify("plasmid: :DnaRevComp [{start}..{stop}|all] (default: the selection)", vim.log.levels.ERROR) end
    if not in_range(rec, s, e) then return end
    se().revcomp(buf, s, e, { force = a.bang })
  end, { nargs = "*", bang = true, desc = "Reverse-complement a range, the selection, or all", complete = function() return { "all" } end })
  cmd("DnaRotate", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    local o = tonumber(a.args)
    if not o and a.args == "" then o = require("plasmid.buffer").cursor_info(vim.api.nvim_get_current_win()) end
    if not o or not in_range(rec, o) then return require("plasmid.util").notify("plasmid: :DnaRotate [{bp}] (default: the cursor's base becomes base 1)", vim.log.levels.ERROR) end
    se().rotate(buf, o)
  end, { nargs = "?", desc = "Move the origin of a circular record" })
  cmd("DnaEdit", function(a)
    local buf, _, rec = cur()
    if not buf then return end
    local s, e
    if a.args:match("^%s*all%s*$") then
      s, e = 1, rec:len()
    else
      s, e = range_args(buf, a.args)
    end
    if not s then
      return require("plasmid.util").notify("plasmid: :DnaEdit [{start}..{stop}|all] (default: the selection)", vim.log.levels.ERROR)
    end
    if not in_range(rec, s, e) then return end
    require("plasmid.scratch").open(buf, s, e)
  end, { nargs = "*", desc = "Edit a range, the selection, or all residues as plain text", complete = function() return { "all" } end })
  cmd("DnaEnzymes", function(a)
    local buf = cur()
    if not buf then return end
    local E = require("plasmid.enzymes")
    local arg = a.fargs[1]
    if not arg then return E.toggle(buf) end
    if arg == "pick" then return E.pick(buf) end
    if arg == "off" or arg == "unique" or arg == "dual" or arg == "all" then return E.set(buf, arg) end
    return E.set(buf, a.fargs)
  end, {
    nargs = "*", desc = "Restriction sites: unique | dual | all | off | pick | Enzyme...",
    complete = function(lead)
      local out = {}
      for _, w in ipairs({ "unique", "dual", "all", "off", "pick" }) do
        if w:find(lead, 1, true) == 1 then out[#out + 1] = w end
      end
      for _, e in ipairs(require("plasmid.bio.restriction").enzymes()) do
        if e.name:lower():find(lead:lower(), 1, true) == 1 then out[#out + 1] = e.name end
      end
      return out
    end,
  })
  cmd("DnaEnzymeImport", function(a)
    local rb = require("plasmid.rebase")
    local R = require("plasmid.bio.restriction")
    local util = require("plasmid.util")
    if a.fargs[1] == "clear" then
      rb.clear(); R.reset()
      return util.notify(("plasmid: imported enzymes removed; using the bundled %d"):format(#R.enzymes()), vim.log.levels.INFO)
    end
    if #a.fargs == 0 then
      local list = rb.load()
      if list then
        return util.notify(("plasmid: using %d imported REBASE enzymes (%s); :DnaEnzymeImport clear reverts"):format(#list,
          vim.fn.fnamemodify(rb.store_path and rb.store_path() or "", ":~")), vim.log.levels.INFO)
      end
      return util.notify(("plasmid: using the bundled %d enzymes; import REBASE with :DnaEnzymeImport emboss_e.NNN [emboss_r.NNN] [emboss_s.NNN]"):format(#R.enzymes()),
        vim.log.levels.INFO)
    end
    local paths = vim.tbl_map(function(p) return vim.fn.expand(p) end, a.fargs)
    local list, sum = rb.import(paths[1], paths[2], paths[3])
    if not list then return util.notify("plasmid: " .. tostring(sum), vim.log.levels.ERROR) end
    R.reset()
    local sk = sum.skipped
    util.notify(("plasmid: imported %d enzymes (%d commercial); skipped %d two-cut, %d without a defined cut, %d bad site, %d duplicate, %d malformed"):format(
      sum.count, sum.commercial, sk.two_cuts or 0, sk.no_cut or 0, sk.bad_site or 0, sk.duplicate or 0, sk.malformed or 0), vim.log.levels.INFO)
    for _, b in ipairs(vim.api.nvim_list_bufs()) do require("plasmid.enzymes").refresh(b) end
  end, { nargs = "*", complete = "file", desc = "Import REBASE EMBOSS files (emboss_e [emboss_r [emboss_s]]) | clear" })
  cmd("DnaFetch", function(a)
    local F = require("plasmid.fetch")
    if a.args == "cancel" then
      local n = F.cancel_all()
      return require("plasmid.util").notify(("plasmid: cancelled %d request%s"):format(n, n == 1 and "" or "s"))
    end
    if a.args == "" then
      return require("plasmid.util").notify("plasmid: :DnaFetch {accession}[:from-to][:-]  (e.g. L09137, NC_000913.3:190-2800)",
        vim.log.levels.WARN)
    end
    F.fetch(a.args, { force = a.bang, refresh = a.bang })
  end, { nargs = "*", bang = true, desc = "Open an NCBI record by accession [:from-to][:-]; ! refreshes and lifts the size limit; cancel" })
  cmd("DnaSearch", function(a)
    require("plasmid.fetch").search(a.args, a.bang)
  end, { nargs = "+", bang = true, desc = "Search NCBI nucleotide (! protein) and open a hit" })
  cmd("DnaGene", function(a)
    require("plasmid.fetch").gene(a.fargs, { force = a.bang })
  end, { nargs = "+", bang = true, desc = "Fetch a gene's genomic region: {symbol} [organism] [flank=N] [plus]" })
  cmd("DnaFetchCache", function(a)
    local nc = require("plasmid.ncbi")
    local util = require("plasmid.util")
    if a.args == "clear" then
      local n = nc.cache_clear()
      return util.notify(("plasmid: removed %d cached record%s"):format(n or 0, n == 1 and "" or "s"))
    end
    local list = nc.cache_list()
    local total = 0
    for _, e in ipairs(list) do total = total + (e.size or 0) end
    util.notify(("plasmid: %d cached record%s, %.1f MB in %s"):format(#list, #list == 1 and "" or "s", total / 1e6,
      nc.config().cache_dir))
  end, { nargs = "?", complete = function() return { "clear" } end, desc = "Show or clear the NCBI record cache" })
  cmd("DnaPrimerAdd", function(a)
    local buf = cur()
    if not buf then return end
    local pr = require("plasmid.primers")
    local first = a.fargs[1]
    if a.range > 0 or first == "fwd" or first == "rev" or not first then
      -- from the selection (visual range or current selection)
      local strand = (first == "rev") and "rev" or "fwd"
      local name = (first == "fwd" or first == "rev") and a.fargs[2] or first
      local sel
      if a.range > 0 then
        local sel_ = require("plasmid.select")
        local s, e = sel_.visual_range(buf)
        if s then sel = sel_.make(require("plasmid.buffer").record(require("plasmid.buffer").get(buf)), s, e) end
      end
      return pr.add_selection(buf, strand, name, sel)
    end
    pr.add_sequence(buf, first, a.fargs[2])
  end, { nargs = "*", range = true, desc = "Add a primer: {sequence} [name], or [fwd|rev] [name] from the selection" })
  cmd("DnaPrimerDesign", function(a)
    local buf = cur()
    if not buf then return end
    require("plasmid.primerdesign").run(buf, a.fargs)
  end, { nargs = "*", desc = "Design a primer pair for the selection: [name=] [tm=] [len=18-30] [slack=8] [dtm=3] [fwd=Enzyme|SEQ] [rev=...] [into=part] [enzyme=BsaI] [homology=25] [overhangs=AATG/GCTT]" })
  cmd("DnaTm", function(a)
    require("plasmid.primerinfo").command(cur(), a.fargs)
  end, { nargs = "*", desc = "Tm (and Ta for two): {seq} [{seq2}] [polymerase=q5] [dmso=%] [conc=nM] [salt=mM]; no seq: selection" })
  cmd("DnaPolymerase", function(a)
    require("plasmid.primerinfo").choose(a.fargs[1])
  end, { nargs = "?", desc = "Polymerase for Tm/Ta (picker without argument)",
    complete = function(lead) return require("plasmid.primerinfo").complete_polymerase(lead) end })
  cmd("DnaPCR", function(a)
    local buf = cur()
    if not buf then return end
    require("plasmid.pcr").run(buf, a.fargs[1], a.fargs[2])
  end, { nargs = "*", desc = "Simulate PCR: two primer names or sequences (default: the record's primers)" })
  cmd("DnaDigest", function(a)
    local buf = cur()
    if not buf then return end
    require("plasmid.digest").run(buf, a.fargs)
  end, { nargs = "*", desc = "Digest with enzymes (default: the named set shown); fragment table and gel",
    complete = function(lead)
      local out = {}
      for _, e in ipairs(require("plasmid.bio.restriction").enzymes()) do
        if e.name:lower():find(lead:lower(), 1, true) == 1 then out[#out + 1] = e.name end
      end
      return out
    end })
  cmd("DnaGel", function(a)
    local buf = cur()
    if not buf then return end
    require("plasmid.digest").gel(buf, a.args)
  end, { nargs = "+", desc = "Gel with one lane per digest: :DnaGel EcoRI | EcoRI HindIII | uncut [ladder=100bp]" })
  -- Cloning simulation ----------------------------------------------------
  -- These run from a sequence buffer or from a digest table, so they take the
  -- current buffer as it is instead of insisting on a record.
  local function clone_parts_complete(lead)
    local out = { "selection" }
    local C = require("plasmid.clone")
    for _, p in ipairs(C.staged) do out[#out + 1] = p.name end
    for b, st in pairs(require("plasmid.buffer").state) do
      if vim.api.nvim_buf_is_valid(b) then out[#out + 1] = require("plasmid.buffer").record(st).name end
    end
    return vim.tbl_filter(function(w) return w:lower():find(lead:lower(), 1, true) == 1 end, out)
  end
  cmd("DnaCloneAdd", function(a)
    require("plasmid.clone").add(vim.api.nvim_get_current_buf(), a.args ~= "" and a.args or nil)
  end, { nargs = "?", desc = "Stage a cloning part: the digest fragment under the cursor, the selection, or this record" })
  cmd("DnaCloneList", function()
    require("plasmid.clone").list(vim.api.nvim_get_current_buf())
  end, { desc = "Show the staged cloning parts (name, length, ends)" })
  cmd("DnaCloneClear", function() require("plasmid.clone").clear() end, { desc = "Empty the staged part list" })
  cmd("DnaClone", function(a)
    require("plasmid.clone").run(vim.api.nvim_get_current_buf(), nil, a.fargs)
  end, { nargs = "*", complete = clone_parts_complete,
    desc = "Assemble parts, inferring the method (Type IIS sites → Golden Gate, terminal homology → Gibson, else ligation)" })
  cmd("DnaLigate", function(a)
    require("plasmid.clone").run(vim.api.nvim_get_current_buf(), "ligation", a.fargs)
  end, { nargs = "*", complete = clone_parts_complete, desc = "Ligate parts by their ends (default: the staged parts or this digest)" })
  cmd("DnaGibson", function(a)
    require("plasmid.clone").run(vim.api.nvim_get_current_buf(), "gibson", a.fargs)
  end, { nargs = "*", complete = clone_parts_complete, desc = "Gibson/HiFi assembly by terminal homology (cloning.min_homology)" })
  cmd("DnaGoldenGate", function(a)
    local args, enzyme = vim.deepcopy(a.fargs), nil
    local first = args[1]
    if first and not require("plasmid.clone").resolve(vim.api.nvim_get_current_buf(), first)
      and require("plasmid.bio.restriction").get(first) then
      enzyme = table.remove(args, 1)
    end
    require("plasmid.clone").run(vim.api.nvim_get_current_buf(), "goldengate", args, { enzyme = enzyme })
  end, { nargs = "*", desc = "Golden Gate assembly: :DnaGoldenGate [enzyme] [parts...] (default enzyme: cloning.enzyme)",
    complete = function(lead)
      local out = clone_parts_complete(lead)
      for _, n in ipairs(require("plasmid.bio.assembly").TYPE_IIS) do
        if n:lower():find(lead:lower(), 1, true) == 1 then out[#out + 1] = n end
      end
      return out
    end })
  cmd("DnaFind", function(a)
    local buf = cur()
    if not buf then return end
    local su = require("plasmid.searchui")
    if a.fargs[1] == nil or a.fargs[1] == "off" then return su.clear(buf) end
    su.run(buf, a.fargs[1], tonumber(a.fargs[2]) or 0, a.bang and -1 or 1)
  end, { nargs = "*", bang = true, desc = "Find a motif (IUPAC; aa:MOTIF for amino acids) [mismatches]; off clears; ! searches backward" })
  cmd("DnaEditReport", function()
    local buf = cur()
    if buf then se().show_report(buf) end
  end, { desc = "Features affected by the last sequence edit" })

  cmd("DnaGroup", function(a)
    local buf = cur()
    if not buf then return end
    local buffer = require("plasmid.buffer")
    local st = buffer.get(buf)
    local v = a.args
    if v == "" then
      -- toggle between continuous and the configured grouping (10 if that is continuous too)
      local conf = tonumber(config.options.linear.group) or 10
      if buffer.group_of(st) == 0 then return buffer.set_group(buf, conf ~= 0 and conf or 10) end
      return buffer.set_group(buf, 0)
    end
    local n = tonumber(v)
    if not n or n < 0 or n ~= math.floor(n) then
      return require("plasmid.util").notify("plasmid: :DnaGroup [n]  (0 = continuous; no argument toggles)", vim.log.levels.ERROR)
    end
    buffer.set_group(buf, n)
  end, { nargs = "?", desc = "Residues per group in the linear view (0 = continuous with a ruler); no argument toggles",
    complete = function() return { "0", "10" } end })
  cmd("DnaTracks", function(a)
    local buf = cur()
    if not buf then return end
    local v = a.args
    if v == "" then return buffer.set_track_limit(buf, nil) end
    if v ~= "all" and not tonumber(v) then
      return require("plasmid.util").notify("plasmid: :DnaTracks {n|all}  (no argument resets to config)", vim.log.levels.ERROR)
    end
    buffer.set_track_limit(buf, v == "all" and "all" or math.max(1, math.floor(tonumber(v))))
  end, { nargs = "?", complete = function() return { "all", "2", "4", "8" } end,
    desc = "Feature rows per line (linear) / rings (map)" })

  cmd("DnaFeatures", function() if cur() then ui.pick_feature() end end, { desc = "Feature picker" })
  cmd("DnaFeatureInfo", function() if cur() then ui.info_under_cursor() end end, { desc = "Feature under cursor" })

  cmd("DnaRecord", function(a)
    local buf, st = cur()
    if not buf then return end
    local function switch(i)
      st.current = i
      st.selection = nil
      require("plasmid.orfs").refresh(buf)
      buffer.render(buf, { keep_bp = 1 })
    end
    if a.args ~= "" then
      local i = tonumber(a.args)
      if not i or not st.records[i] then
        return require("plasmid.util").notify(("plasmid: record must be 1..%d"):format(#st.records), vim.log.levels.ERROR)
      end
      return switch(i)
    end
    local idx = {}
    for i = 1, #st.records do idx[i] = i end
    vim.ui.select(idx, {
      prompt = "Record",
      format_item = function(i)
        local r = st.records[i]
        return ("%d  %s  %d bp  %s"):format(i, r.name, r:len(), r.description:sub(1, 50))
      end,
    }, function(i) if i then switch(i) end end)
  end, { nargs = "?", desc = "Select record in a multi-record file" })

  cmd("DnaTopology", function(a)
    local buf, st, rec = cur()
    if not buf then return end
    if a.args ~= "linear" and a.args ~= "circular" then
      return require("plasmid.util").notify("plasmid: :DnaTopology linear|circular", vim.log.levels.ERROR)
    end
    if rec:is_protein() then
      return require("plasmid.util").notify("plasmid: protein sequences have no topology", vim.log.levels.WARN)
    end
    if rec.topology == a.args then return end
    local tx = require("plasmid.model.transaction")
    tx.commit(rec, { op = "topology", value = a.args }, "topology " .. a.args)
    if rec.format == "fasta" then
      -- FASTA has no topology field: persist it as an NCBI-style header modifier
      local d = rec.description:gsub("%s*%[topology=%a+%]", ""):gsub("%s*%(circular%)", "")
      if a.args == "circular" then d = (d .. " [topology=circular]"):gsub("^%s+", "") end
      rec.description = d
      rec.raw.header_dirty = true
    end
    buffer.refresh_modified(buf)
    if a.args == "linear" and st.view == "circular" then st.view = "linear" end
    buffer.render(buf)
  end, { nargs = 1, complete = function() return { "linear", "circular" } end, desc = "Set topology" })

  cmd("DnaInfo", function()
    local buf, st, rec = cur()
    if not buf then return end
    local alphabet = require("plasmid.model.alphabet")
    local counts = {}
    for _, f in ipairs(rec.features) do counts[f.type] = (counts[f.type] or 0) + 1 end
    local types = vim.tbl_keys(counts)
    table.sort(types)
    local lines = {
      ("%s  (%s)"):format(rec.name, rec.id or ""),
      rec.description,
      rec:is_protein() and ("%d aa  ·  protein"):format(rec:len())
        or ("%d %s  ·  %s  ·  %s  ·  GC %.1f%%"):format(rec:len(), rec:unit(), rec.molecule_type, rec.topology,
          100 * alphabet.gc(rec.seq:tostring())),
      ("file: %s  (%s, record %d/%d)"):format(vim.fn.fnamemodify(st.path, ":~:."), st.format, st.current, #st.records),
      "",
    }
    for _, t in ipairs(types) do lines[#lines + 1] = ("  %-18s %d"):format(t, counts[t]) end
    vim.lsp.util.open_floating_preview(lines, "", { border = "rounded", focus_id = "plasmid_info" })
  end, { desc = "Record summary" })

  cmd("DnaRaw", function()
    local buf, st = cur()
    if not buf then return end
    local fh = io.open(st.path, "rb")
    if not fh then return require("plasmid.util").notify("plasmid: cannot read " .. st.path, vim.log.levels.ERROR) end
    local text = fh:read("*a")
    fh:close()
    if st.format == "snapgene" then text = require("plasmid.io.snapgene").describe(text) end
    vim.cmd("vsplit")
    local raw = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(raw, 0, -1, false, vim.split(text, "\r?\n"))
    vim.bo[raw].modifiable = false
    vim.bo[raw].filetype = ({ genbank = "genbank", fasta = "fasta", snapgene = "xml" })[st.format]
    vim.api.nvim_buf_set_name(raw, "plasmid-raw://" .. st.path)
    require("plasmid.util").case_insensitive_search(raw)
    vim.api.nvim_win_set_buf(0, raw)
  end, { desc = "Show the file's raw text (read-only)" })

  cmd("DnaCalibrate", function()
    -- draws a circle with the current aspect so the user can judge roundness
    local Canvas = require("plasmid.view.canvas")
    local aspect = config.options.circular.aspect
    local w, h = 40, 20
    local cv = Canvas.new(w, h)
    local r = 4 * h / 2 - 3
    for k = 0, 720 do
      local th = 2 * math.pi * k / 720
      cv:dot(w + r * math.sin(th) / aspect, 2 * h - r * math.cos(th))
    end
    local lines = cv:render()
    table.insert(lines, 1, ("circular.aspect = %.2f   (circle looks tall → raise it; wide → lower it)"):format(aspect))
    vim.lsp.util.open_floating_preview(lines, "", { border = "rounded", focus_id = "plasmid_cal" })
  end, { desc = "Check circle aspect ratio for this font" })
end

local function register_autocmds()
  local buffer = require("plasmid.buffer")
  group = vim.api.nvim_create_augroup("plasmid", { clear = true })
  local o = config.options
  if o.autocmds then
    vim.api.nvim_create_autocmd("BufReadCmd", {
      group = group, pattern = o.patterns,
      callback = function(ev)
        buffer.open(ev.buf, vim.fn.fnamemodify(ev.match, ":p"))
      end,
    })
    vim.api.nvim_create_autocmd({ "BufWriteCmd", "FileWriteCmd" }, {
      group = group, pattern = o.patterns,
      callback = function(ev)
        if buffer.get(ev.buf) then
          buffer.write(ev.buf, vim.fn.fnamemodify(ev.match, ":p"))
        else
          -- not ours (e.g. a plain buffer saved under a sequence-file name)
          vim.cmd(("noautocmd write! %s"):format(vim.fn.fnameescape(ev.match)))
          vim.bo[ev.buf].modified = false
        end
      end,
    })
  end
  vim.api.nvim_create_autocmd({ "BufWipeout", "BufDelete" }, {
    group = group,
    callback = function(ev) buffer.detach(ev.buf) end,
  })
  vim.api.nvim_create_autocmd({ "WinResized", "VimResized" }, {
    group = group,
    callback = function()
      local seen = {}
      for _, win in ipairs(vim.api.nvim_list_wins()) do
        local b = vim.api.nvim_win_get_buf(win)
        local st = buffer.get(b)
        if st and st.layout and not seen[b] then
          seen[b] = true
          local w = vim.api.nvim_win_get_width(win)
          local h = vim.fn.winheight(win)
          local lay = st.layout
          if math.abs((lay.width or 0) - w) > 2 or (lay.view == "circular" and lay.height ~= h) then
            vim.api.nvim_win_call(win, function() buffer.render(b) end)
          end
        end
      end
    end,
  })
  vim.api.nvim_create_autocmd("BufWinEnter", {
    group = group,
    callback = function(ev)
      local st = buffer.get(ev.buf)
      if st and config.options.winbar then
        vim.wo[0][0].winbar = "%{%v:lua.require'plasmid.buffer'.winbar()%}"
        vim.wo[0][0].wrap = false
        vim.wo[0][0].list = false
        vim.wo[0][0].spell = false
        vim.wo[0][0].cursorline = true
      end
    end,
  })
  vim.api.nvim_create_autocmd("CursorMoved", {
    group = group,
    callback = function(ev)
      if buffer.get(ev.buf) then require("plasmid.codon").update(ev.buf) end
    end,
  })
  vim.api.nvim_create_autocmd("ColorScheme", {
    group = group,
    callback = function() require("plasmid.view.highlights").setup() end,
  })
end

function M.setup(opts)
  config.setup(opts)
  require("plasmid.view.highlights").setup()
  require("plasmid.buffer").setup_decorations()
  register_commands()
  register_autocmds()
  initialized = true
end

--- Called from plugin/ at startup; a later setup() call reconfigures.
function M._autoload()
  if not initialized then M.setup(vim.g.plasmid or {}) end
end

return M
