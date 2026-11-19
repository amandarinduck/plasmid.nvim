local M = {}

M.defaults = {
  -- Take over reading/writing of sequence files (BufReadCmd/BufWriteCmd).
  autocmds = true,
  patterns = {
    "*.gb", "*.gbk", "*.gbff", "*.genbank", "*.ape",
    "*.fa", "*.fasta", "*.fna", "*.ffn", "*.fas", "*.fsa",
    "*.faa", "*.pep", "*.gp", "*.gpff",   -- protein FASTA / GenPept
    "*.dna", "*.prot", "*.rna",           -- SnapGene
  },
  -- Use feature colours stored in the file (SnapGene segments, ApE/Benchling
  -- /ApEinfo_*color) when they are legible; otherwise colour by feature type.
  use_file_colors = true,
  -- Open reading frames shown by `go` / :DnaOrfs
  orfs = {
    min_aa = 75,              -- SnapGene's default minimum
    starts = "atg",           -- "atg" or "alternative" (all start codons of the genetic code)
    table = nil,              -- genetic code; nil = record's /transl_table, else 11 (circular) / 1
    hide_annotated = true,    -- hide ORFs ending at the stop codon of an annotated CDS
  },
  -- Plain yanks (no register given) also go to the system clipboard:
  -- true = "+" register, "auto" = follow Neovim's 'clipboard' option, false = never.
  -- Without a clipboard provider the copy is sent to the terminal as OSC 52.
  yank_clipboard = true,
  -- Sequence editing
  -- restriction enzymes (:DnaEnzymes)
  enzymes = {
    default_set = "unique",   -- what gz / :DnaEnzymes shows: "unique" | "dual" | "all"
    min_site = 6,             -- sets skip enzymes with fewer specific (non-N) site letters
    lanes = 4,                -- label rows above each sequence line (more are summarised as +n)
    commercial_only = true,   -- with an imported REBASE set: sets skip enzymes without suppliers
  },
  -- primers (:DnaPrimerAdd, :DnaTm, :DnaPCR)
  primers = {
    polymerase = "q5",        -- Tm/Ta model (:DnaPolymerase lists them; NEB Tm Calculator rules)
    conc_nM = nil,            -- primer concentration override (default: the polymerase's)
    salt_mM = nil,            -- monovalent salt override (default: the polymerase buffer's)
    dmso = 0,                 -- % DMSO (lowers Ta by 0.6 °C per %)
    show = true,              -- draw primer_bind features as primer bases in the linear view
    lanes = 3,                -- primer rows above / below each sequence line
    binding = { anchor = 10, max_mismatches = 3, resume = 5, min_anneal = 10 },
    design = { len = { 18, 30 }, tm = { 57, 68 }, max_dtm = 3, gc = { 0.35, 0.65 }, max_run = 4, pairs = 8 },
  },
  -- NCBI access (:DnaFetch, :DnaSearch, :DnaGene, gF). Needs curl.
  ncbi = {
    api_key = nil,            -- NCBI API key (or $NCBI_API_KEY): 10 requests/s instead of 3
    email = nil,              -- sent to NCBI only if you set it (their usage policy asks for one)
    max_length = 5000000,     -- larger records need :DnaFetch!
    cache = true,             -- keep fetched records in stdpath("cache")/plasmid/ncbi
    cache_days = 30,
  },
  fetch = {
    dir = nil,                -- where fetched records are proposed to be saved (default: cwd)
    search_results = 20,      -- hits listed by :DnaSearch
    gene_flank = 0,           -- bases added on each side by :DnaGene (flank=N overrides)
    organism = nil,           -- default organism for :DnaGene, e.g. "human"
  },
  -- simulated gels (:DnaDigest, :DnaGel)
  gel = {
    ladder = "auto",          -- "auto" | "1kb" | "100bp"
    height = 24,              -- rows from the wells to the bottom
  },
  -- cloning simulation (:DnaClone, :DnaLigate, :DnaGoldenGate, :DnaGibson)
  cloning = {
    min_homology = 20,        -- shortest terminal overlap Gibson/HiFi will join
    max_homology = 200,       -- longest overlap looked for (a longer match is reported at this length)
    enzyme = "BsaI",          -- default Type IIS enzyme for :DnaGoldenGate and for method inference
    annotate_parts = true,    -- mark each part in the product as a misc_feature
    annotate_junctions = true,-- mark each junction (overhang or overlap) in the product
    max_solutions = 50,       -- stop the search after this many assemblies (the rest are not offered)
  },
  edit = {
    -- features cut by a deletion: "truncate" (just report), "flag" (also add
    -- a /note to the feature), "refuse" (cancel the edit; ! forces it)
    cut_policy = "truncate",
    -- rewrite /translation of CDS features whose bases changed (removed when
    -- the coding sequence now has an internal stop codon)
    update_translations = true,
  },
  -- Feature detection (:DnaDetect) against imported databases + your library
  detect = {
    min_identity = 0.95,
    min_coverage = 0.95,      -- of the reference, for a full-length hit
    fragments = true,         -- also propose partial matches...
    min_fragment_coverage = 0.5,
    min_fragment_length = 50,
    nested_same_type = false, -- hide hits inside a better hit of the same type
  },
  featuredb_dir = nil,        -- nil = stdpath('data')/plasmid/featuredb
  -- Personal feature library (GenBank, one record per feature); nil = stdpath('data')/plasmid/library.gb
  library_path = nil,
  -- "linear" | "circular" | "auto" (auto = circular for circular records)
  default_view = "linear",
  -- Qualifiers consulted, in order, for a feature's display label.
  label_qualifiers = { "label", "gene", "product", "locus_tag", "standard_name", "allele", "note" },
  -- Feature types never drawn.
  hidden_types = { source = true },
  -- Which features win space in crowded regions. Score = type score
  -- (+ labelled bonus when the feature has a real name; + pinned bonus for
  -- features you added in this session). Ties: longer first.
  priority = {
    types = {
      CDS = 100, gene = 95, mRNA = 90, ncRNA = 85, tRNA = 85, rRNA = 85,
      rep_origin = 80, oriT = 80, promoter = 80, terminator = 75,
      RBS = 70, polyA_signal = 70, enhancer = 65, LTR = 65, regulatory = 60,
      protein_bind = 60, primer_bind = 55, sig_peptide = 55, mat_peptide = 55,
      exon = 50, misc_RNA = 50, ORF = 45, mobile_element = 45, intron = 40,
      misc_feature = 30, misc_binding = 25, repeat_region = 20, rep_element = 20,
      variation = 10, misc_difference = 5, old_sequence = 5, conflict = 5, unsure = 5,
    },
    default = 40,
    labelled_bonus = 10,
    pinned_bonus = 1000,
  },
  -- Hide a `gene` whose location is identical to a CDS/mRNA/tRNA/rRNA.
  dedupe_genes = true,
  -- Show the record's winbar (name, length, cursor bp, feature under cursor).
  winbar = true,
  linear = {
    line_width = "auto",     -- bases per line, or "auto" to fit the window (multiple of 10)
    show_complement = true,
    show_translation = true,  -- per-CDS translation, in the CDS's own frame and code
    color_bases = false,
    case = "auto",            -- "auto" (upper for GenBank, preserve for FASTA) | "upper" | "preserve"
    max_tracks = 4,           -- feature rows per line; the rest are summarised as "+N more"
    group = 10,               -- residues per group; 0 = continuous with a ruler row (codons never split by gaps)
  },
  circular = {
    -- Physical width/height of one braille dot. Braille cells are 2x4 dots and
    -- terminal cells are ~1:2, so dots are ~square. Raise this if circles look
    -- tall, lower it if they look wide (:DnaCalibrate helps).
    aspect = 1.0,
    max_tracks = 6,
    label_width = 22,
    show_ticks = true,
    -- additional types hidden on the map only (annotation history, not biology)
    hidden_types = { old_sequence = true, misc_difference = true, conflict = true, variation = true, unsure = true },
  },
  -- Mouse in a sequence buffer ('mouse' must be enabled in Neovim).
  mouse = {
    enabled = true,
    drag_select = true,          -- press and drag to select a region
    double_click_feature = true, -- double click selects the feature under the pointer
    primer_pair = true,          -- Ctrl-click two primers to select the product they make
  },
  keymaps = {
    toggle_view = "<Tab>",
    activate = "<CR>",
    info = "K",
    next_feature = "]f",
    prev_feature = "[f",
    features = "gO",
    feature_outer = "af",     -- text object: whole feature span (yaf, vaf)
    feature_inner = "if",     -- text object: feature's own residues, spliced, 5'->3' (yif, vif)
    yank_bases = "y",         -- visual: yank bases under the selection (top strand)
    yank_bases_rc = "Y",      -- visual: yank reverse complement (bottom strand 5'->3')
    yank_bases_pick = "gy",   -- visual: choose the strand and direction to yank
    yank_selection = "gy",    -- yank the current selection
    yank_selection_rc = "gY",
    clear_selection = "<Esc>",
    toggle_orfs = "go",       -- show/hide open reading frames
    add_feature = "ga",       -- add a feature over the selection (normal) or visual range
    edit_feature = "gE",      -- edit the feature under the cursor
    follow_link = "gF",       -- open a record linked from the feature (protein_id, GeneID, UniProt)
    delete_feature = "",      -- unmapped by default; :DnaFeatureDelete
    undo = "u",
    redo = "<C-r>",
    -- sequence editing (in the plasmid view)
    insert_before = "i",      -- insert residues before the cursor base (prompt)
    insert_after = "a",       -- ... after it
    delete_base = "x",        -- delete [count] bases; visual: delete the selection
    delete_visual = "d",      -- visual: delete the selection (daf/dif delete a feature)
    replace_base = "r",       -- [count]r{letter}: replace bases
    change_visual = "c",      -- visual: replace the selection (caf: a feature)
    paste_after = "p",        -- insert register contents after the cursor (visual: replace)
    paste_before = "P",
    revcomp_visual = "R",     -- visual: reverse-complement the selection in place
    -- motif search
    search_forward = "/",     -- prompt for a motif (IUPAC; "aa:" for amino acids; "MOTIF 1" = 1 mismatch)
    search_backward = "?",
    search_next = "n",
    search_prev = "N",
    select_hit = "gn",        -- select the current hit (then gy, ga, c, ...)
    add_primer = "gp",        -- visual: make a primer from the selection (forward or reverse)
    clone_add = "gA",         -- stage this record (or the selection) as a cloning part (:DnaCloneAdd)
    toggle_enzymes = "gz",    -- show/hide restriction sites (enzymes.default_set)
    pick_enzymes = "gZ",      -- choose enzymes from a list with site counts
    expand_tracks = "+",      -- show every feature track
    collapse_tracks = "-",    -- back to the configured track limit
  },
  -- Foreground colours per feature type (hex). Missing types use `palette`.
  colors = {
    CDS = "#f2b134", gene = "#86c06c", mRNA = "#7fa7d9", exon = "#9ab8e0",
    promoter = "#4fbf75", terminator = "#e05d5d", rep_origin = "#e3d24a", oriT = "#e3d24a",
    primer_bind = "#4fb3d9", protein_bind = "#b07cd8", RBS = "#de8f4e", regulatory = "#5fbfa0",
    polyA_signal = "#d070a0", enhancer = "#8fd16a", LTR = "#c09060", misc_feature = "#a8a8a8",
    sig_peptide = "#e090c0", mat_peptide = "#e0b070", repeat_region = "#8c8c8c",
    mobile_element = "#c0a080", ORF = "#8fa3bf", tRNA = "#70c0c0", rRNA = "#50a0a0", ncRNA = "#80b0b0",
    misc_RNA = "#80b0b0",
  },
  palette = { "#e6a23c", "#67c23a", "#409eff", "#f56c6c", "#b37feb", "#36cfc9", "#ff9c6e", "#bae637" },
}

M.options = vim.deepcopy(M.defaults)

function M.setup(opts)
  M.options = vim.tbl_deep_extend("force", vim.deepcopy(M.defaults), opts or {})
  return M.options
end

return M
