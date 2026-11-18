--- The fixed set of NCBI requests recorded into tests/fixtures/ncbi_http/ by
--- tests/tools/ncbi_record.lua, replayed by tests/spec/ncbi_spec.lua and
--- re-run against NCBI with PLASMID_LIVE=1.
---
--- kind = "fetch"   spec, opts, cfg (ncbi option overrides)
---        "summary" db, ids, cfg
---        "search"  db, term, opts
---        "gene"    symbol, organism
---        "eutil"   util, params (raw E-utility call, response kept as is)
local C = {}

C.list = {
  { name = "L09137", kind = "fetch", spec = "L09137", note = "pUC19c; overlaps tests/fixtures/ncbi/L09137.gb" },
  { name = "J01749", kind = "fetch", spec = "J01749", note = "pBR322; overlaps tests/fixtures/ncbi/J01749.gb" },
  { name = "NC_001422", kind = "fetch", spec = "NC_001422", note = "phiX174; overlaps tests/fixtures/ncbi/NC_001422.gb" },
  { name = "NM_000518.5", kind = "fetch", spec = "NM_000518.5", note = "HBB mRNA" },
  { name = "thr_plus", kind = "fetch", spec = "NC_000913.3:190-2800", note = "E. coli thr operon region, plus strand" },
  { name = "thr_minus", kind = "fetch", spec = "NC_000913.3:2800-190", note = "same region, reverse complement" },
  { name = "NP_000509.1", kind = "fetch", spec = "NP_000509.1", note = "HBB protein (GenPept)" },
  { name = "P69905", kind = "fetch", spec = "P69905", note = "UniProt/Swiss-Prot accession (HBA_HUMAN) in db protein" },
  { name = "con_LN625376", kind = "fetch", spec = "LN625376.1", note = "CON scaffold record (1 contig); gbwithparts" },
  { name = "gene_thrA_region", kind = "fetch", spec = "NC_000913.3:337-2799",
    note = "the region :DnaGene thrA fetches (gene coordinates from gene_thrA)" },
  { name = "unknown", kind = "fetch", spec = "XX999999", expect_err = true, note = "unknown accession" },
  { name = "big_refused", kind = "fetch", spec = "NC_000913.3", cfg = { max_length = 1000000 }, expect_err = true,
    note = "E. coli genome (4.6 Mb) refused by the size guard: only the esummary is made" },
  { name = "gi_number", kind = "fetch", spec = "20141090", note = "GI of L09137.2; efetch shared with L09137" },
  { name = "summary_post", kind = "summary", db = "nuccore", ids = { "L09137", "J01749", "NC_001422" },
    cfg = { post_ids = 2 }, note = "esummary sent as POST (id list longer than post_ids)" },
  { name = "search_puc19", kind = "search", db = "nuccore", term = "pUC19[title] AND complete sequence[title]" },
  { name = "search_protein", kind = "search", db = "protein",
    term = "hemoglobin subunit beta[title] AND Homo sapiens[orgn] AND refseq[filter]", opts = { retmax = 5 } },
  { name = "search_nohit", kind = "search", db = "nuccore", term = "zzqxqzzq[title]", note = "no hits: warnings" },
  { name = "gene_BRCA1", kind = "gene", symbol = "BRCA1", organism = "Homo sapiens", note = "minus-strand gene" },
  { name = "gene_thrA", kind = "gene", symbol = "thrA", organism = "Escherichia coli str. K-12 substr. MG1655",
    note = "plus-strand gene, inside the thr_plus range" },
  { name = "raw_con_gb", kind = "eutil", util = "efetch",
    params = { { "db", "nuccore" }, { "id", "LN625376.1" }, { "rettype", "gb" }, { "retmode", "text" } },
    note = "the same CON record with rettype=gb: CONTIG line, no sequence" },
  { name = "raw_efetch_unknown", kind = "eutil", util = "efetch",
    params = { { "db", "nuccore" }, { "id", "XX999999" }, { "rettype", "gbwithparts" }, { "retmode", "text" } },
    note = "efetch error format (HTTP 200, letter-spaced text)" },
  { name = "raw_gene_672", kind = "eutil", util = "efetch",
    params = { { "db", "gene" }, { "id", "672" }, { "retmode", "text" } },
    note = "NCBI's own BRCA1 gene record (Annotation line) to check the derived coordinates" },
  { name = "raw_gene_945803", kind = "eutil", util = "efetch",
    params = { { "db", "gene" }, { "id", "945803" }, { "retmode", "text" } },
    note = "NCBI's own thrA gene record" },
}

C.by_name = {}
for _, c in ipairs(C.list) do C.by_name[c.name] = c end

--- Run one case; returns the handle. Config overrides are applied for the
--- duration of the call setup (the client reads options when a call starts).
function C.run(ncbi, case, cb, extra_cfg)
  local cfgmod = require("plasmid.config")
  local saved = cfgmod.options.ncbi
  cfgmod.options.ncbi = vim.tbl_extend("force", saved or {}, extra_cfg or {}, case.cfg or {})
  local h
  if case.kind == "fetch" then
    h = ncbi.fetch(case.spec, case.opts or {}, cb)
  elseif case.kind == "summary" then
    h = ncbi.summary(case.db, case.ids, cb)
  elseif case.kind == "search" then
    h = ncbi.search(case.db, case.term, case.opts or {}, cb)
  elseif case.kind == "gene" then
    h = ncbi.gene(case.symbol, case.organism, case.opts or {}, cb)
  elseif case.kind == "eutil" then
    h = ncbi.eutil(case.util, case.params, {}, cb)
  end
  cfgmod.options.ncbi = saved
  return h
end

return C
