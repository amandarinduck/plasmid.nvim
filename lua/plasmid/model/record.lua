-- Record and Feature model. Pure Lua.
local Store = require("plasmid.model.store")
local location = require("plasmid.model.location")

local M = {}

---------------------------------------------------------------------------
-- Feature
---------------------------------------------------------------------------

local Feature = {}
Feature.__index = Feature
M.Feature = Feature

--- @param t { type: string, location_text: string, qualifiers?: table, raw_lines?: string[] }
function Feature.new(t)
  local f = setmetatable({
    id = t.id,
    type = t.type,
    location_text = t.location_text,
    qualifiers = t.qualifiers or {}, -- ordered { {key=, value=, quoted=bool}, ... }
    raw_lines = t.raw_lines,         -- verbatim source lines, nil once modified
    color = t.color,                 -- display colour from the file (SnapGene), if any
    no_arrow = t.no_arrow,           -- non-directional feature (or protein record)
    origin = t.origin,               -- "primer" for SnapGene primer sites
    raw_xml = t.raw_xml,             -- verbatim SnapGene <Feature> element, nil once modified
    recent_id = t.recent_id,         -- SnapGene feature id
    primer_seq = t.primer_seq,       -- full primer 5'->3' (SnapGene primers), incl. tails
    snapgene_site = t.snapgene_site, -- SnapGene's annealing info for this primer site, if read
  }, Feature)
  f:_parse_location()
  return f
end

function Feature:_parse_location()
  local node, err = location.parse(self.location_text)
  self.loc = node
  self.loc_error = err
  if node then
    self.parts = location.flatten(node)
    self.strand = location.strand(self.parts)
  else
    self.parts = {}
    self.strand = 1
  end
end

--- First value of a qualifier, or nil.
function Feature:get(key)
  for _, q in ipairs(self.qualifiers) do
    if q.key == key then return q.value end
  end
  return nil -- explicit: callers pass this straight to tonumber()
end

--- All values for a qualifier key.
function Feature:get_all(key)
  local out = {}
  for _, q in ipairs(self.qualifiers) do
    if q.key == key then out[#out + 1] = q.value end
  end
  return out
end

function Feature:mark_dirty()
  self.raw_lines = nil
  self.raw_xml = nil
end

local DEFAULT_LABEL_KEYS = { "label", "gene", "product", "locus_tag", "standard_name", "allele", "note" }

--- Display label, following SnapGene/ApE convention (/label first).
function Feature:label(keys)
  for _, k in ipairs(keys or DEFAULT_LABEL_KEYS) do
    local v = self:get(k)
    if v and v ~= true and v ~= "" then
      v = tostring(v):gsub("%s+", " ")
      return v
    end
  end
  return self.type
end

--- Colour stored in the file: SnapGene segment colour, or ApE/Benchling
--- /ApEinfo_fwdcolor and /ApEinfo_revcolor qualifiers.
function Feature:file_color()
  if self.color then return self.color end
  local q = self.strand == -1 and self:get("ApEinfo_revcolor") or self:get("ApEinfo_fwdcolor")
  if type(q) == "string" and q:match("^#%x%x%x%x%x%x$") then return q end
  return nil
end

function Feature:length()
  return location.length(self.parts)
end

--- 5' and 3' end positions on the forward axis (nil for remote-only).
function Feature:bounds()
  local lo, hi
  for _, p in ipairs(self.parts) do
    if not p.ref then
      lo = lo and math.min(lo, p.start) or p.start
      hi = hi and math.max(hi, p.stop) or p.stop
    end
  end
  return lo, hi
end

---------------------------------------------------------------------------
-- Record
---------------------------------------------------------------------------

local Record = {}
Record.__index = Record
M.Record = Record

--- @param t table fields: name, id, description, molecule_type, topology,
---   division, date, seq (string), features, format, raw (format-specific)
function Record.new(t)
  local r = setmetatable({
    name = t.name or t.id or "unnamed",
    id = t.id or t.name,
    description = t.description or "",
    molecule_type = t.molecule_type or "DNA",
    alphabet = t.alphabet or "dna",   -- "dna" | "rna" | "protein"
    topology = t.topology or "linear",
    division = t.division,
    date = t.date,
    seq = Store.new(t.seq or ""),
    features = {},
    format = t.format,
    raw = t.raw or {},
    revision = 0,
  }, Record)
  for _, f in ipairs(t.features or {}) do r:add_feature(f) end
  return r
end

function Record:is_protein()
  return self.alphabet == "protein"
end

--- Unit for positions: "bp" (DNA), "nt" (RNA), "aa" (protein).
function Record:unit()
  return ({ dna = "bp", rna = "nt", protein = "aa" })[self.alphabet] or "bp"
end

function Record:len()
  return self.seq:len()
end

function Record:is_circular()
  return self.topology == "circular" and self.alphabet ~= "protein"
end

function Record:add_feature(f)
  if getmetatable(f) ~= Feature then f = Feature.new(f) end
  self._next_id = (self._next_id or 0) + 1
  f.id = self._next_id
  if self.alphabet == "protein" then f.no_arrow = true end
  self.features[#self.features + 1] = f
  return f
end

function Record:feature_by_id(id)
  if id and id < 0 then return self.overlay_by_id and self.overlay_by_id[id] or nil end
  for _, f in ipairs(self.features) do
    if f.id == id then return f end
  end
end

--- Extents of a feature on this record (handles origin wrap).
function Record:extents(f)
  return location.extents(f.parts, self:len(), self:is_circular())
end

--- Residues of a feature in biological orientation.
function Record:extract(f)
  return location.extract(f.parts, self.seq)
end

--- Change topology. Marks the header dirty so writers regenerate it.
function Record:set_topology(topo)
  assert(topo == "linear" or topo == "circular", "topology must be linear|circular")
  assert(not (topo == "circular" and self.alphabet == "protein"), "a protein cannot be circular")
  if self.topology ~= topo then
    self.topology = topo
    self.raw.locus_dirty = true
    self.revision = self.revision + 1
  end
end

return M
