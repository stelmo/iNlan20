function getmet(model, metid)
    for met in model["metabolites"]
        if met["id"] == metid
            return met
        end
    end
    return Dict()
end

function getgene(model, geneid)
    for gene in model["genes"]
        if gene["id"] == geneid
            return gene
        end
    end
    return Dict()
end

function splitformula(str)
    indstarts = findall(isuppercase, str)
    indends = [indstarts[2:end] .- 1; length(str)]
    f = Dict()
    for (a, b) in zip(indstarts, indends)
        p = str[a:b]
        ind = findfirst(isdigit, p)
        if isnothing(ind)
            f[p[1:end]] = 1.0
        else
            f[p[1:ind-1]] = parse(Float64, p[ind:end])
        end
    end
    return f
end

function getrxn(model, rxnid)
    for rxn in model["reactions"]
        if rxn["id"] == rxnid
            return rxn
        end
    end
    return Dict()
end

function remrxns(model, rxnstoremove)
    rxnlist = []
    for rxn in model["reactions"]
        if rxn["id"] in rxnstoremove
            continue
        end
        push!(rxnlist, rxn)
    end
    model["reactions"] = rxnlist
end

function remmets(model, metstoremove)
    metlist = []
    for met in model["metabolites"]
        if met["id"] in metstoremove
            continue
        end
        push!(metlist, met)
    end
    model["metabolites"] = metlist
end

function remo2rxns(model)
    rxnlist = []
    for rxn in model["reactions"]
        if "o2_c" in keys(rxn["metabolites"])
            continue
        end
        push!(rxnlist, rxn)
    end
    model["reactions"] = rxnlist
end


function curategenes!(model)
    genes = String[] # get all genes associated with reactions
    for k=1:length(model["reactions"])
        if model["reactions"][k]["gene_reaction_rule"] != ""
            append!(genes, split(model["reactions"][k]["gene_reaction_rule"], r"(\s)(or|and)(\s)"))
        end
    end

    newgeneset = [] # add all genes with reaction associations
    for pid in genes
        gene = Dict()
        gene["id"] = pid
        gene["name"] = "Protein ID: $(pid)"
        gene["notes"] = Dict()
        gene["annotation"] = Dict("sbo" => ["SBO:0000243"], "jgi-protein-id" => pid)
        push!(newgeneset, gene)
    end
    model["genes"] = newgeneset


    # remove repeated genes in grr  - won't be repeats in the manually added ones
    for k=1:length(model["reactions"])
        grrs = unique(split(model["reactions"][k]["gene_reaction_rule"], " or "))
        model["reactions"][k]["gene_reaction_rule"] = join(grrs, " or ")
    end

    # remove duplicate genes
    alreadyadded = String[]
    genelist = []
    for k=1:length(model["genes"])
        id = model["genes"][k]["id"]
        if id in alreadyadded
            continue
        end
        push!(genelist, model["genes"][k])
        push!(alreadyadded, id)
    end
    model["genes"] = genelist

end

"""
fix gene reaction rule
"""
function fixgrr!(model, rxnid, newgpr, iseccode=false)
    if iseccode
        for k=1:length(model["reactions"])
            if rxnid in get(model["reactions"][k]["annotation"], "ec-code", [])
                # println(rxnid)
                model["reactions"][k]["gene_reaction_rule"] = newgpr
            end
        end
    else
        for k=1:length(model["reactions"])
            if model["reactions"][k]["id"] == rxnid
                model["reactions"][k]["gene_reaction_rule"] = newgpr
            end
        end
    end
end


"""
Ensure the annotations are in the correct format for compatibility
"""
function fixannotations!(bigg, bigg_ec=Dict())
    for k in keys(bigg)
        tmpdict = Dict{String, Any}()

        # Assign SBO terms
        if haskey(bigg[k], "metabolites")
            tmpdict["bigg.reaction"] = [k]
            if occursin("EX_", k)
                tmpdict["sbo"] = "SBO:0000627"
            elseif length(k) > 3 && occursin("t", k[end-3:end])
                tmpdict["sbo"] = "SBO:0000185"
            else
                tmpdict["sbo"] = "SBO:0000176"
            end
        else
            tmpdict["bigg.metabolite"] = [k]
            tmpdict["sbo"] = "SBO:0000247"
        end

        for anno in bigg[k]["annotation"]

            if string(anno[1]) == "BioCyc" # nothing fancy here
                if haskey(tmpdict, "biocyc")
                    push!(tmpdict["biocyc"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["biocyc"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("MetaNetX", string(anno[1]))
                if occursin("Equation", string(anno[1]))
                    if haskey(tmpdict, "metanetx.reaction")
                        push!(tmpdict["metanetx.reaction"], split(string(anno[2]),"/")[end])
                    else
                        tmpdict["metanetx.reaction"] = [split(string(anno[2]),"/")[end]]
                    end
                else
                    if haskey(tmpdict, "metanetx.chemical")
                        push!(tmpdict["metanetx.chemical"], split(string(anno[2]),"/")[end])
                    else
                        tmpdict["metanetx.chemical"] = [split(string(anno[2]),"/")[end]]
                    end
                end
            end

            if "RHEA" == string(anno[1])
                if haskey(tmpdict, "rhea")
                    push!(tmpdict["rhea"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["rhea"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if "InChI Key" == string(anno[1])
                if haskey(tmpdict, "inchikey")
                    push!(tmpdict["inchikey"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["inchikey"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("SEED Compound", string(anno[1]))
                if haskey(tmpdict, "seed.compound")
                    push!(tmpdict["seed.compound"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["seed.compound"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("SEED Reaction", string(anno[1]))
                if haskey(tmpdict, "seed.reaction")
                    push!(tmpdict["seed.reaction"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["seed.reaction"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("CHEBI", string(anno[1]))
                if haskey(tmpdict, "chebi")
                    push!(tmpdict["chebi"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["chebi"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("Reactome", string(anno[1]))
                if haskey(tmpdict, "reactome")
                    push!(tmpdict["reactome"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["reactome"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("KEGG Compound", string(anno[1]))
                if haskey(tmpdict, "kegg.compound")
                    push!(tmpdict["kegg.compound"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["kegg.compound"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if occursin("KEGG Reaction", string(anno[1]))
                if haskey(tmpdict, "kegg.reaction")
                    push!(tmpdict["kegg.reaction"], split(string(anno[2]),"/")[end])
                else
                    tmpdict["kegg.reaction"] = [split(string(anno[2]),"/")[end]]
                end
            end

            if haskey(bigg_ec, k) # use own ec assignments preferentially
                tmpdict["ec-code"] = bigg_ec[k]
            else
                if occursin("EC Number", string(anno[1]))
                    if haskey(tmpdict, "ec-code")
                        push!(tmpdict["ec-code"], split(string(anno[2]),"/")[end])
                    else
                        tmpdict["ec-code"] = [split(string(anno[2]),"/")[end]]
                    end
                end
            end
        end

        bigg[k]["annotation"] = tmpdict
    end
end

"""
newmodel = newmodel()

Creates a blank model.
"""
function newmodel(modelname)
    model = Dict()
    model["genes"] = []
    model["version"] = "0"
    model["reactions"] = Array{Dict, 1}()
    model["compartments"] = Dict("c"=>"cytosol","e"=>"extracellular space","h"=>"hydrogenosome")
    model["id"] = modelname
    model["metabolites"] = Array{Dict, 1}()
    return model
end


"""
universal_X = getuniX(db, X)

Returns a dictionary where the ID points to the data for X = {reactions, metabolites}
"""
function getuniX(db, X)
    uniX = Dict()
    for x in db[X]
        uniX[x["id"]] = x
    end
    return uniX
end

"""
mets = getmets(rxn)

Return a list of all the metabolites used in a reaction
"""
function getmets(rxndict)
    mets = String[]
    for met in keys(rxndict["metabolites"])
        push!(mets, met)
    end
    return mets
end

"""
addrxns!(model, data, dbref)

Add reactions to model.
NB: Make all reactions reversible.
"""
function addrxns!(model, data, dbref)
    added = Array{Dict, 1}()
    for dpoint in unique(data)
        tmprxn = Dict() # don't add annotations
        tmprxn["id"] = dbref[dpoint]["id"]
        tmprxn["name"] = dbref[dpoint]["name"]
        tmprxn["notes"] = dbref[dpoint]["notes"]
        tmprxn["metabolites"] = dbref[dpoint]["metabolites"]
        tmprxn["gene_reaction_rule"] = dbref[dpoint]["gene_reaction_rule"]
        tmprxn["annotation"] = dbref[dpoint]["annotation"]
        # Make reversible
        tmprxn["lower_bound"] = -1000.0
        tmprxn["upper_bound"] = 1000.0
        push!(added, tmprxn)
    end
    model["reactions"] = added
end

"""
addmets!(model, data, dbref)

Add metabolites to the model.
"""
function addmets!(model, data, dbref)
    added = Array{Dict, 1}()
    for dpoint in unique(data)
        tmpmet = Dict() # don't add annotations
        tmpmet["id"] = dbref[dpoint]["id"]
        tmpmet["name"] = dbref[dpoint]["name"]
        tmpmet["notes"] = dbref[dpoint]["notes"]
        tmpmet["annotation"] = dbref[dpoint]["annotation"]
        tmpmet["charge"] = get(dbref[dpoint], "charge", 0)
        tmpmet["formula"] = get(dbref[dpoint], "formula", "X")

        # default metabolites go to the cellular compartment
        if endswith(tmpmet["id"], "_e")
            tmpmet["compartment"] = "e"
            if haskey(dbref, dpoint[1:end-1]*"c")
                tmpmet["charge"] = dbref[dpoint[1:end-1]*"c"]["charge"]
                tmpmet["formula"] = dbref[dpoint[1:end-1]*"c"]["formula"]
            end
        else
            tmpmet["compartment"] = "c"
        end
        push!(added, tmpmet)
    end
    model["metabolites"] = added
end

"""
addtemprxn!(metabolite_id)

Adds a reversible reaction to the model
that takes in metabolite_id from nothing and
supplies or removes it from the system

Note, the metabolites should already be present
in the model.
"""
function addtemprxn!(model, met)
    rxnid = string("TMP_", met)
    rxn = Dict()
    rxn["name"] = string("Temporary ", met, " exchange")
    rxn["metabolites"] = Dict(met=>-1)
    rxn["lower_bound"] = -1000.0
    rxn["upper_bound"] = 1000.0
    rxn["id"] = rxnid
    rxn["annotation"] = Dict()
    rxn["gene_reaction_rule"] = ""
    rxn["notes"] = Dict()
    push!(model["reactions"], rxn)
end

"""
changebound!(model, rxn, lb, ub)
"""
function changebound!(model, rxnid, lb, ub)
    found = false
    for rxn in model["reactions"]
        if rxn["id"] == rxnid
            rxn["lower_bound"] = lb
            rxn["upper_bound"] = ub
            found = true
            break
        end
    end
    !found && @info("Reaction $(rxnid) not found.")
end

"""
Adds or updates the BOF
"""
function addupbof!(model, metdict)
    bof = "Biomass"
    found = false
    for rxn in model["reactions"]
        if rxn["id"] == bof
            rxn["metabolites"] = metdict
            found = true
            break
        end
    end
    if !found
        rxn = Dict()
        rxn["name"] = string("Biomass objective function")
        rxn["metabolites"] = metdict
        rxn["lower_bound"] = 0.0
        rxn["upper_bound"] = 1000.0
        rxn["id"] = bof
        rxn["annotation"] = Dict("sbo" => "SBO:0000629")
        rxn["gene_reaction_rule"] = ""
        rxn["notes"] = Dict("description" => "Biomass objective function reaction")
        push!(model["reactions"], rxn)
    end
end

"""
rmtmprxns(model)

Removes all temporary reactions from the model
"""
function rmtmprxns!(model)
    model["reactions"] = filter(x->!startswith(x["id"],"TMP"), model["reactions"])
end

"""
addsink(model, met)

Add a sink reaction (takes metabolite and converts it to nothing)
"""
function addsink!(model, met)
    rxnid = string("R_SINK_", met)
    rxn = Dict()
    rxn["name"] = string(met, " sink")
    rxn["metabolites"] = Dict(met=>-1)
    rxn["lower_bound"] = 0
    rxn["upper_bound"] = 1000.0
    rxn["id"] = rxnid
    rxn["annotation"] = Dict("sbo" => "SBO:0000632")
    rxn["gene_reaction_rule"] = ""
    rxn["notes"] = Dict("description"=>"Sink reaction")
    push!(model["reactions"], rxn)
end

"""
addsupply(model, met)

Add a demand reaction (takes metabolite and converts it to nothing)
"""
function addsupply!(model, met)
    rxnid = string("R_DM_", met)
    rxn = Dict()
    rxn["name"] = string(met, " supply")
    rxn["metabolites"] = Dict(met=>-1)
    rxn["lower_bound"] = -1000.0
    rxn["upper_bound"] = 0.0
    rxn["id"] = rxnid
    rxn["annotation"] = Dict("sbo" => "SBO:0000628")
    rxn["gene_reaction_rule"] = ""
    rxn["notes"] = Dict("description" => "Supply/demand reaction")
    push!(model["reactions"], rxn)
end

"""
addrxn(model, rxndict)

Add a reaction only if it isn't already in the model
"""
function addrxn!(model, rxndict)
    found = false
    for rxn in model["reactions"]
        if rxn["id"] == rxndict["id"]
            found = true
            break
        end
    end
    !found && push!(model["reactions"], rxndict)
end

"""
addmet(model, metdict)

Add a metabolite only if it isn't already in the model
"""
function addmet!(model, metdict)
    found = false
    for met in model["metabolites"]
        if met["id"] == metdict["id"]
            found = true
            break
        end
    end
    !found && push!(model["metabolites"], metdict)
end

"""
"""
function createmetabolite(id, name, compartment="c", charge=0, formula="", anno=Dict())
    met = Dict()
    met["id"] = id
    met["name"] = name
    met["notes"] = Dict()
    met["charge"] = charge
    met["formula"] = formula
    met["annotation"] = Dict()
    met["compartment"] = compartment
    return met
end

function createmetabolite(hid, unimets)
    met = copy(unimets[hid[1:end-1]*"c"])
    met["id"] = hid
    met["compartment"] = "h"
    return met
end

function createreaction(id, name, metdict; lb=-1000, ub=1000, anno=Dict(), notes=Dict(), gpr="")
    rxn = Dict()
    rxn["id"] = id
    rxn["name"] = name
    rxn["metabolites"] = metdict
    rxn["lower_bound"] = lb
    rxn["upper_bound"] = ub
    rxn["annotation"] = Dict()
    rxn["gene_reaction_rule"] = gpr
    rxn["notes"] = notes
    return rxn
end

function fixmisc(model)
    # add general metabolite SBO terms
    metids = ["fdxox_h","fdxrd_h","rq_ox_h","rq_red_h","xylan_e","hicit_c","s2ab_c","glglutl2abut_c","opth_c","cellulose_e","hemicellulose_e","acoa_c","cer_18_p_c","sphmye_c","ppchol_c","spppchol_c","gluccera_c","lacsylcera_c","galsylcera_c","digalsylcera_c"]
    for metid in metids
        met = getmet(model, metid)
        if !haskey(met, "annotation")
            met["annotation"] = Dict()
        end
        met["annotation"]["sbo"] = "SBO:0000247"
    end

    # add metabolic reaction sbo term
    rxnids = ["HYDh3","FUMh","r_xylanase","CMPL2h","SUCCOAh","SUCOAACTh","PFLh","r_dextrinase","HYDhbi","PFOh","r_HAH","MEhx","MASh","r_gapd","HYDhfe","MEhy","NGAM","r_R04019","r_R01100","r_R02543","r_R01498","r_hemicellulase","r_R06522","r_R01500","r_lacd1","r_lacd2","r_02541","HYDh2","H2ht","r_glutabutl","r_2abut2akgt","r_HICITD","r_cellulase","r_MetSmethtrans","r_glglutglyl"]
    for rxnid in rxnids
        # println(rxnid)
        rxn = getrxn(model, rxnid)
        if !haskey(rxn, "annotation")
            rxn["annotation"] = Dict()
        end
        rxn["annotation"]["sbo"] = "SBO:0000176"
    end

    trxnids = ["CYSabc","ARB_Dabc","H2Oht","RAFFabc","AKGMAL","SUCCabc","GLNabc","MALTabc","ARGabc","MALht","SUCCht2","ARBabc","ATPShydr","ACht2","PIht","SUCCht1","Kt1","CO2ht","PYRht","PPIabc","CMPL1h","CELLBabc","SERabc","ADPATPht","RMNabc","SUCRabc","FORht2","NMNTP","NAt"]
    for trxnid in trxnids
        trxn = getrxn(model, trxnid)
        if !haskey(trxn, "annotation")
            trxn["annotation"] = Dict()
        end
        trxn["annotation"]["sbo"] = "SBO:0000185"
    end

    for rxnid in ["HYDh2", "HYDh3"]
        trxn = getrxn(model, rxnid)
        trxn["annotation"]["kegg.reaction"] = "K18332"

    end

    for metid in ["nh3_c", "metox_c"]
        met = getmet(model, metid)
        delete!(met["annotation"], "inchikey")
    end

    remrxns(model, ["ANNAT", "NADDP","ASPO5_1", "METSOXR2","METSOXR1","PAPSR"])

    # remove metabolites
    mets2rem = ["o2_c","n4abutn_c","trypta_c","48dhoxquin_c","fe3_c","oaa_h","akg_h","glu__L_h","asp__L_h", "pap_c"]
    remmets(model, mets2rem)

end
