function genes_uniprot(db)
    gu = Dict()
    for (k, v) in db
        for vv in split(v[2], " ")
            gu[vv] = k
        end
    end
    return gu
end

function reversebidir(bidir)
    rbidir = Dict()
    for (k, v) in bidir
        rbidir[v] = k
    end
    return rbidir
end

function getcazygenes(loc)
    cellulases = [5,6,7,8,9,12,44,45,48,51,74,124]
    hemicellulases = [5,8,10,11,16,26,30,43,44,51,62,98,141]
    cazygenes = Dict()
    open(loc) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            cprts = [split(x, "_")[1] for x in split(prts[2],"-")]
            cazygenes[prts[1]] = cprts
        end
    end
    cls = String[]
    hls = String[]
    for (k, cs) in cazygenes
        for cid in cellulases
            dom = string("GH", cid)
            if dom in cs
                push!(cls, k)
                break
            end
        end

        for cid in hemicellulases
            dom = string("GH", cid)
            if dom in cs
                push!(hls, k)
                break
            end
        end
    end
    return cls, hls
end

function get_ec_genes(loc, jgipos, bidirpos, onlybidir=false)

    ec_genes = Dict()
    open(loc) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            id = prts[1]
            if onlybidir
                jgi_ec = []
            else
                jgi_ec = [prts[jgipos]]
            end
            bidir_ec = split(prts[bidirpos], "; ")
            ecs = [x for x in unique([jgi_ec; bidir_ec]) if x != ""]
            if !isempty(ecs)
                for ec in ecs
                    if haskey(ec_genes, ec)
                        push!(ec_genes[ec], id)
                    else
                        ec_genes[ec] = [id]
                    end
                end
            end
        end
    end

    return ec_genes

end

function kegg_maps(loc)
    kmaps = Dict()

    allfile = open(loc) do io
        read(io, String)
    end
    lns = split(allfile, "\\n")

    for ln in lns
        prts = split(ln, "\\t")
        if length(prts) > 1
            kmaps[prts[2]] = prts[1]
        end
    end

    return kmaps
end

function readecs(loc)
    ecs = String[]
    open(loc) do io
        for ln in eachline(io)
            push!(ecs, ln)
        end
    end
    return ecs
end

function getsubsystems(modelfileloc, modelsheets)
    rxn_subs = Dict()

    for sheet in modelsheets
        db = readxlsheet(modelfileloc, sheet)
        for k=2:size(db, 1)
            rxn_subs[db[k, 2]] = db[k, 4]
        end
    end

    return rxn_subs
end

"""
Reads the TCDB file
"""
function readtransporterfile(tfloc)
    transporters = Dict()
    open(tfloc) do f
        for ln in eachline(f)
            prts = split(ln, "\t")
            pid = prts[1]
            tc = prts[4]
            transporters[pid] = tc
        end
    end
    return transporters
end

"""
dict = parse_ipr_pfam(floc, start_and_stop_pos=false)

Returns a dictionary mapping a protein id => [(domain name1, starting position1), (domain name2, starting position2)...]
using ONLY PFam annotations.
"""
function parse_ipr_pfam(floc)
    pid_annos = Dict()
    firstline = true
    open(floc) do f
        for line in readlines(f)
            firstline && (firstline=false; continue)
            parts = split(line, "\t")
            pid = parts[1]
            db = parts[4]
            dom = parts[5] # 6 is the description
            if db == "HMMPfam"
                numhits = parse(Int64, parts[7])
                for k=1:numhits
                    pos = parse(Int64, split(parts[8], ",")[k])
                    posend = parse(Int64, split(parts[9], ",")[k])
                    if haskey(pid_annos, pid)
                        push!(pid_annos[pid], dom)
                    else
                        pid_annos[pid] = [dom]
                    end
                end
            end
        end
    end
    return pid_annos
end


"""
ecannos = getECannos(ecloc)

Return the EC numbers associated with protein ids. Return only a single EC number per protein ID
by returning the EC with the least number of dashes (the most complete one)
"""
function getECannos(ecloc)
    ecdict = Dict()
    open(ecloc) do f
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            proteinid = prts[1]
            ecnum = prts[2]
            ecdict[proteinid] = push!(get(ecdict, proteinid, String[]), ecnum)
        end
    end
    ecs = Dict()
    for (k, v) in ecdict
         uv = unique(v)
        if length(uv) > 1
            numdash = 10
            for uvx in uv
                @info string(ecloc, " ", k, " had multiple assignments.")
                n = count(x->x=="-", split(uvx, "."))
                if n < numdash
                    numdash = n
                    ecs[k] = uvx
                end
            end
        else
            ecs[k] = uv[1]
        end
    end
    return ecs
end

"""
ids = matchids(ploc, pos1=4, pos2=3)

Matches transcript IDs to protein IDs
"""
function matchids(ploc, pos1=4, pos2=3)
    ids = Dict()
    reader = FASTA.Reader(open(ploc, "r"))
    for record in reader
        idfull = split(BioSequences.FASTA.identifier(record), "|")
        id1 = idfull[pos1]
        id2 = idfull[pos2]
        haskey(ids, id1) && @warn("ERROR: duplicate id.")
        ids[id1] = id2

    end
    close(reader)
    return ids
end

"""
tevid = getTranscriptomicEvidence(tloc)

Return the protein IDs (data from mkdbfiles) that had a transcript associated with it.
"""
function getTranscriptomicEvidence(tloc, evalcutoff = 1e-60, qcovcutoff = 90.0)
    tevid = Dict()
    open(tloc) do f
        for ln in eachline(f)
            prts = split(ln, "\t")
            transcriptid = prts[1] # from the transcriptome not JGI
            proteinid = prts[2]
            eval = parse(Float64, prts[3])
            qcov = parse(Float64, prts[4])
            if eval < evalcutoff && qcov > qcovcutoff
                if haskey(tevid, proteinid)
                    push!(tevid[proteinid], transcriptid)
                else
                    tevid[proteinid] = [transcriptid]
                end
            end
        end
    end
    return tevid
end


"""
pid_to_uid = matchbidir(fungus_to_uniprot_loc, uniprot_to_fungus_loc, evalcutoff, ecposition)

Returns a dictionary of the bidirectional hits meeting evalcutoff.
"""
function matchbidir(fungus_to_uniprot_loc, uniprot_to_fungus_loc, evalcutoff=1e-20, idpos=3)

    fungus_uni = readbidirblast(fungus_to_uniprot_loc, evalcutoff)
    uni_fungus = readbidirblast(uniprot_to_fungus_loc, evalcutoff)
    bidir = Dict()
    for (k, v) in fungus_uni
        if get(uni_fungus, v, "X") == k
            if idpos == 3 # normal genome input
                bidir[k] = split(v, "|")[2] # only return the entry id
            else
                bidir[k] = v # only return the entry id
            end
        end
    end

    return bidir
end

"""
bidirdict = readbidirblast(fileloc, evalcutoff)

Return a dictionary mapping the query to the db (q => db) in a dict, selecting the highest coverage (first) db hit
and checking that the evalue is below the cutoff
"""
function readbidirblast(floc, evalcutoff, format=-1, evpid=false; kmod=false)
    bidir = Dict()
    open(floc) do f
        for ln in eachline(f)
            prts = split(ln, "\t")
            if kmod
                query = split(prts[1], "|")[2]
            else
                query = prts[1]
            end
            if format == -1
                db = prts[2]
            else
                db = split(prts[2], "|")[format]
            end
            eval = parse(Float64, prts[3])
            if eval < evalcutoff && !haskey(bidir, query) # highest cov is shown first
                if evpid
                    bidir[query] = [db, eval, parse(Float64, prts[4])]
                else
                    bidir[query] = db
                end
            end
        end
    end
    return bidir
end

"""
entry_descr_mapping = readuniprot(floc)

Returns the entry mapping from the uniprot description file
"""
function readuniprot(floc, retids=[1,4,5])
    uni = Dict()
    open(floc) do f
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            entry = prts[1]
            uni[entry] = [prts[x] for x in retids]
        end
    end
    return uni
end


"""
entry_descr_mapping = readtcdb(floc)

Returns the entry mapping from the uniprot description file
"""
function readtcdb(floc)
    tcdb = Dict()
    open(floc) do f
        for ln in eachline(f)
            prts = split(ln, "\t")
            entry = prts[1]
            tcid = prts[2]
            tcdesc = prts[3]
            tcdb[entry] = [tcid, tcdesc]
        end
    end
    return tcdb
end

"""
ecs = getecsmaster(masterdict)

Returns a list of unique EC numbers associated with a master annotation dictionary
"""
function getecsmaster(mdict)
    ecdict = Dict()
    for nit in keys(mdict)
        jgiec = mdict[nit][2]
        tid = mdict[nit][3]
        haskey(ecdict, jgiec) ? nothing : (ecdict[jgiec]=[])
        push!(ecdict[jgiec], tid) # no transcript => "" entry (but has gene associated!NB)
        uniecs = split(mdict[nit][5], "; ")
        for ec in uniecs
            haskey(ecdict, ec) ? nothing : (ecdict[ec]=[])
            push!(ecdict[ec], tid)
        end
    end
    for k in keys(ecdict)
        ecdict[k] = unique(ecdict[k])
    end
    return ecdict
end

function getecsmasterS3(mdict)
    ecdict = Dict()
    for tid in keys(mdict)
        uniecs = split(mdict[tid][3], "; ")
        for ec in uniecs
            haskey(ecdict, ec) ? nothing : (ecdict[ec]=[])
            push!(ecdict[ec], tid)
        end
    end
    for k in keys(ecdict)
        ecdict[k] = unique(ecdict[k])
    end
    return ecdict
end

function gettcdbmaster(loc)
    master = Dict()
    open(loc) do f
        headings = ["Protein ID","Transcript ID", "TCID"]
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            pid = prts[1]
            tid = prts[2]
            tcid = prts[3]
            master[tcid] = tid # "" => has gene but no transcript
        end
    end
    return master

end

"""
allbiggrxns = readExcelModel(excelloc)

Return a list of bigg reactions IDs from all the sheets in loc
"""
function readExcelModel(loc, blocks)
    rxns = String[]
    ecs = Dict()
    for block in blocks
        data = readxlsheet(loc, block)
        for (i, t) in enumerate(data[2:end, 2])
            if !isna(t)
                push!(rxns, convert(String, t))
                if haskey(ecs, data[i+1, 1])
                    push!(ecs[data[i+1, 1]], convert(String, t))
                else
                    ecs[data[i+1, 1]] = [convert(String, t)]
                end
            end
        end
    end

    return ecs, unique(filter(x->x!="", rxns))
end

"""
ec_genes = get_gene_ecs(masterannotationdict)

Returns the gene -> ec list association
"""
function get_gene_ecs(mat)
    ec_genes = Dict()
    for (tid, entry) in mat
        jgiec = [entry[2]]
        uecs = split(entry[5], "; ")
        ecs = filter(x->x!="", [jgiec; uecs])
        for ec in ecs
            if haskey(ec_genes, ec)
                push!(ec_genes[ec], tid)
            else
                ec_genes[ec] = [tid]
            end
        end
    end
    return ec_genes
end

"""
eclist = readkegglink(filelocation)

Reads the link file from Kegg and returns the associated list of ECs.
"""
function readkegglink(floc)
    ecs = String[]
    open(floc) do f
        for ln in eachline(f)
            prts = split(ln, "\t")
            ec = split(prts[2], "ec:")[2]
            push!(ecs, ec)
        end
    end
    return unique(ecs)
end

"""
records_dict = readall(file_location, is_aminoacid, id_loc=3)

Read in all the sequence data and push to array records. id_loc = the position of the id for this entry
"""
function readall(ploc, isaa::Bool=false, idloc=3, delim="|")
    records = Dict()
    reader = FASTA.Reader(open(ploc, "r"))
    if isaa
        for record in reader
            idfull = FASTA.identifier(record)
            idname = split(idfull, delim)[idloc]
            records[idname] = FASTA.sequence(LongAminoAcidSeq, record)
        end
    else
        for record in reader
            idfull = FASTA.identifier(record)
            idname = split(idfull, delim)[idloc]
            records[idname] = FASTA.sequence(LongDNASeq, record)
        end
    end
    close(reader)
    return records
end

"""
writeall(floc, seqs)

Write the seqs to floc
"""
function writeall(floc, seqs)
    w = FASTA.Writer(open(floc, "w"))
    for (k, v) in seqs
        rec = FASTA.Record(k, v)
        write(w, rec)
    end
    close(w)
end

"""
materanno = readmaster(floc)

Reads the master annotation files
"""
function readmastermetabolic(floc)
    master = Dict()
    open(floc) do f
        headings = ["Transcript ID", "Protein ID", "JGI EC", "Transcriptomic Match", "Proteomics Coverage", "Bidir EC", "Bidir Subunit"]
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            tid = prts[1]
            pid = prts[2]
            jgiec = prts[3]
            transcript = prts[4]
            pcov = prts[5]
            bidirecs = prts[6]
            bidirsubunit = prts[7]
            master[tid] = [pid, jgiec, transcript, pcov, bidirecs, bidirsubunit]
        end
    end
    return master
end

function readmastermetabolicS3(floc)
    master = Dict()
    open(floc) do f
        headings = ["Transcript ID", "Protein ID", "JGI EC", "Transcriptomic Match", "Proteomics Coverage", "Bidir EC", "Bidir Subunit"]
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            tid = prts[1]
            uid = prts[2]
            ecs = prts[3]
            master[tid] = [tid, uid, ecs]
        end
    end
    return master
end

function readmastertransporters(floc)
    master = Dict()
    open(floc) do f
        headings = ["Protein ID", "Transcriptomic Match", "TCID"]
        firstline = true
        for ln in eachline(f)
            firstline && (firstline=false; continue)
            prts = split(ln, "\t")
            pid = prts[1]
            tid = prts[2]
            tcid = prts[3]
            master[pid] = [tid, tcid]
        end
    end
    return master
end


"""
writemodel(model)

Writes the model dict to a json file.
"""
function savemodel(model, modelloc, suppressinfo=false)
    if !suppressinfo
        @info "Saving model with..."
        @info string(length(model["reactions"])," reactions")
        @info string(length(model["metabolites"])," metabolites")
        @info string(length(model["genes"])," genes")
    end
    open(modelloc, "w") do file
        JSON.print(file, model)
    end
end

"""
loadmodel(modelloc)

Loads the model file
"""
function loadmodel(modelloc)
    model = open(modelloc) do file
        JSON.parse(file)
    end
    return model
end

"""
ecs = readECfile(loc)

Read the file containing all the ECs for a fungus
"""
function readECfile(loc)
    ecs = String[]
    open(loc) do f
        for ln in eachline(f)
            push!(ecs, ln)
        end
    end
    return unique(filter(x->x!="", ecs))
end

"""
Read the expression data - only returns
"""
function readexpressiondata(floc, prefix="CB")
    firstline = true
    cols = String[]
    d = Dict()
    open(floc) do f
        for line in eachline(f)
            if firstline
                cols = split(line, "\t")
                firstline = false
            else
                prts = split(line, "\t")
                id = prts[1]
                dd = Float64[]
                for (tpm, col) in zip(prts[2:end], cols)
                    if startswith(col, prefix)
                        push!(dd, parse(Float64, tpm))
                    end
                end
                d[id] = [mean(dd), std(dd)] # assume 3 replicates
            end
        end
    end
    return d
end

function readexpressiondataS3(floc, cond="M2")
    d = Dict()
    firstline = true
    open(floc) do f
        for line in eachline(f)
            if firstline
                firstline = false
                continue
            else
                line = replace(line,"\"" => "")
                line = replace(line, "\\" => "")
                prts = split(line, ",")
                id = prts[1]
                cd = prts[2]
                startswith(cd, cond) ? nothing : continue
                if haskey(d, id)
                    d[id] += parse(Float64, prts[4])
                else
                    d[id] = parse(Float64, prts[4])
                end
            end
        end
    end
    for k in keys(d)
        d[k] = d[k]/3.0
    end
    return d
end
