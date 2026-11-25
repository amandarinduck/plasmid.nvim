-- Highlight groups. All are defined with `default = true` so colourschemes
-- and user config override them.
local config = require("plasmid.config")

local M = {}

local function cterm(hex)
  -- nearest xterm-256 colour-cube entry, for 'notermguicolors'
  local r, g, b = tonumber(hex:sub(2, 3), 16), tonumber(hex:sub(4, 5), 16), tonumber(hex:sub(6, 7), 16)
  local function q(v) return v < 48 and 0 or (v < 115 and 1 or math.floor((v - 35) / 40)) end
  return 16 + 36 * q(r) + 6 * q(g) + q(b)
end

local function group_name(ftype)
  return "PlasmidFeature" .. ftype:gsub("[^%w]", "_")
end

local defined = {}
local file_groups = {}

local function set(name, spec)
  spec.default = true
  vim.api.nvim_set_hl(0, name, spec)
end

function M.setup()
  defined = {}
  file_groups = {}
  set("PlasmidTitle", { link = "Title" })
  set("PlasmidHeader", { link = "Comment" })
  set("PlasmidRuler", { link = "LineNr" })
  set("PlasmidSequence", { link = "Normal" })
  set("PlasmidComplement", { link = "Comment" })
  set("PlasmidBackbone", { link = "Normal" })
  set("PlasmidTick", { link = "LineNr" })
  set("PlasmidLeader", { link = "NonText" })
  set("PlasmidOverflow", { link = "Comment" })
  set("PlasmidSelection", { link = "Visual" })
  set("PlasmidCodon", { link = "IncSearch" })
  set("PlasmidSearch", { link = "Search" })
  set("PlasmidSearchCurrent", { link = "CurSearch" })
  set("PlasmidSearchInfo", { link = "Search" })
  set("PlasmidEnzyme", { fg = "#d7875f", ctermfg = 173 })
  set("PlasmidEnzymeUnique", { fg = "#ff875f", ctermfg = 209, bold = true })
  set("PlasmidEnzymeCut", { fg = "#ff5f5f", ctermfg = 203, bold = true })
  set("PlasmidEnzymeSite", { underline = true, sp = "#ff875f" })
  set("PlasmidGelLadder", { fg = "#8a8a8a", ctermfg = 245 })
  set("PlasmidGelBand", { fg = "#d0d0d0", ctermfg = 252 })
  set("PlasmidGelBandStrong", { fg = "#ffffff", ctermfg = 231, bold = true })
  set("PlasmidGelLabel", { link = "Number" })
  set("PlasmidGelWell", { fg = "#5f5f5f", ctermfg = 59 })
  set("PlasmidPrimer", { fg = "#4fb3d9", ctermfg = 74, bold = true })
  set("PlasmidPrimerTail", { fg = "#87afaf", ctermfg = 109, italic = true })
  set("PlasmidPrimerMismatch", { fg = "#ffffff", bg = "#d75f5f", ctermfg = 231, ctermbg = 167, bold = true })
  set("PlasmidPrimerArrow", { fg = "#4fb3d9", ctermfg = 74, bold = true })
  set("PlasmidPrimerLabel", { fg = "#4fb3d9", ctermfg = 74 })
  set("PlasmidCodonInfo", { link = "Special" })
  -- alignment view (:DnaAlign, :DnaAlignReads, :DnaCompare)
  set("PlasmidAlignMatch", { link = "Comment" })
  set("PlasmidAlignMismatch", { link = "DiffText" })
  set("PlasmidAlignGap", { link = "DiffDelete" })
  set("PlasmidAlignCoord", { link = "LineNr" })
  set("PlasmidAlignSilent", { link = "Comment" })
  set("PlasmidAlignMissense", { link = "WarningMsg" })
  set("PlasmidAlignNonsense", { link = "ErrorMsg" })
  -- whole-record overview (:DnaOverview). The cells themselves are coloured
  -- with the feature groups above; these are the furniture around them.
  set("PlasmidOverviewCoord", { link = "LineNr" })
  set("PlasmidOverviewAxis", { link = "NonText" })
  set("PlasmidOverviewEmpty", { link = "NonText" })
  set("PlasmidOverviewCount", { link = "Comment" })
  set("PlasmidOverviewLabel", { link = "Comment" })
  set("PlasmidOverviewViewport", { fg = "#ffcc00", ctermfg = 220, bold = true })
  set("PlasmidOverviewCursor", { fg = "#14161b", bg = "#ffcc00", ctermfg = 16, ctermbg = 220, bold = true })
  -- genome atlas (the circular map above circular.atlas.threshold). The two
  -- density rings borrow the gene and mRNA colours so a reader coming from a
  -- plasmid map recognises the strands; GC and skew are a warm/cool pair each.
  set("PlasmidRingBase", { link = "NonText" })
  set("PlasmidRingForward", { fg = "#86c06c", ctermfg = 107 })
  set("PlasmidRingReverse", { fg = "#7fa7d9", ctermfg = 110 })
  set("PlasmidRingGCHigh", { fg = "#e3d24a", ctermfg = 184 })
  set("PlasmidRingGCLow", { fg = "#8a7f3a", ctermfg = 100 })
  set("PlasmidRingSkewPos", { fg = "#4fb3d9", ctermfg = 74 })
  set("PlasmidRingSkewNeg", { fg = "#d87c7c", ctermfg = 174 })
  set("PlasmidRingOri", { fg = "#ff875f", ctermfg = 209, bold = true })
  set("PlasmidLandmark", { fg = "#ffcc00", ctermfg = 220, bold = true })
  set("PlasmidProposed", { fg = "#c8a2e8", ctermfg = 183, italic = true })
  set("PlasmidSelectionArc", { fg = "#ffcc00", ctermfg = 220, bold = true })
  set("PlasmidCenter", { link = "Title" })
  set("PlasmidBaseA", { fg = "#6cc070", ctermfg = 71 })
  set("PlasmidBaseC", { fg = "#5aa0e6", ctermfg = 74 })
  set("PlasmidBaseG", { fg = "#d6b84a", ctermfg = 179 })
  set("PlasmidBaseT", { fg = "#e06c6c", ctermfg = 167 })
  set("PlasmidBaseN", { link = "Comment" })
  -- amino-acid classes (Clustal-like), used when linear.color_bases is on
  set("PlasmidAAHydrophobic", { fg = "#80a0f0", ctermfg = 111 })
  set("PlasmidAAPositive", { fg = "#f01505", ctermfg = 196 })
  set("PlasmidAANegative", { fg = "#c048c0", ctermfg = 170 })
  set("PlasmidAAPolar", { fg = "#15c015", ctermfg = 40 })
  set("PlasmidAACysteine", { fg = "#f08080", ctermfg = 210 })
  set("PlasmidAAGlycine", { fg = "#f09048", ctermfg = 209 })
  set("PlasmidAAProline", { fg = "#c0c000", ctermfg = 184 })
  set("PlasmidAAAromatic", { fg = "#15a4a4", ctermfg = 37 })
  set("PlasmidAAStop", { link = "ErrorMsg" })
  for ftype, hex in pairs(config.options.colors) do
    local name = group_name(ftype)
    set(name, { fg = hex, ctermfg = cterm(hex), bold = ftype == "CDS" })
    defined[ftype] = name
  end
  for i, hex in ipairs(config.options.palette) do
    set("PlasmidPalette" .. i, { fg = hex, ctermfg = cterm(hex) })
  end
end

--- Highlight group for a feature type.
function M.for_type(ftype)
  if defined[ftype] then return defined[ftype] end
  local h = 0
  for i = 1, #ftype do h = (h * 31 + ftype:byte(i)) % 2147483647 end
  return "PlasmidPalette" .. (h % #config.options.palette + 1)
end

--- Highlight group for a feature: its colour from the file (SnapGene/ApE)
--- when enabled and legible on both light and dark backgrounds, else by type.
function M.for_feature(f)
  if config.options.use_file_colors then
    local hex = f.file_color and f:file_color()
    if hex and hex:match("^#%x%x%x%x%x%x$") then
      local r, g, b = tonumber(hex:sub(2, 3), 16), tonumber(hex:sub(4, 5), 16), tonumber(hex:sub(6, 7), 16)
      local lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255
      if lum > 0.12 and lum < 0.92 then
        local name = "PlasmidFile_" .. hex:sub(2):lower()
        if not file_groups[name] then
          vim.api.nvim_set_hl(0, name, { fg = hex, ctermfg = cterm(hex), default = true })
          file_groups[name] = true
        end
        return name
      end
    end
  end
  return M.for_type(f.type)
end

local AA_CLASS = {
  A = "Hydrophobic", I = "Hydrophobic", L = "Hydrophobic", M = "Hydrophobic", F = "Hydrophobic",
  W = "Hydrophobic", V = "Hydrophobic", K = "Positive", R = "Positive", D = "Negative", E = "Negative",
  N = "Polar", Q = "Polar", S = "Polar", T = "Polar", C = "Cysteine", G = "Glycine", P = "Proline",
  H = "Aromatic", Y = "Aromatic", ["*"] = "Stop",
}
M.residue_groups = {}
for k, v in pairs(AA_CLASS) do
  M.residue_groups[k] = "PlasmidAA" .. v
  M.residue_groups[k:lower()] = "PlasmidAA" .. v
end

M.base_groups = {
  A = "PlasmidBaseA", C = "PlasmidBaseC", G = "PlasmidBaseG", T = "PlasmidBaseT", U = "PlasmidBaseT",
  a = "PlasmidBaseA", c = "PlasmidBaseC", g = "PlasmidBaseG", t = "PlasmidBaseT", u = "PlasmidBaseT",
}

return M
