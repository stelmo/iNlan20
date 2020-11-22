#################################################################
# Model of the hydrogenosome
#################################################################

# Metabolites
met = GSM.createmetabolite("pyr_h", unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("coa_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("accoa_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("succ_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("succoa_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("atp_h", unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("adp_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("pi_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("ac_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("h_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("h2_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("co2_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("mal__L_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("nad_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("nadh_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("nadp_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("nadph_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("oaa_h", unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("akg_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("glu__L_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("asp__L_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("for_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("fum_h",  unimets)
GSM.addmet!(model, met)
met = GSM.createmetabolite("h2o_h",  unimets)
GSM.addmet!(model, met)

met = GSM.createmetabolite("fdxox_h", "Oxidized ferredoxin hydrogenosome", "h", 1, "FeS", Dict("kegg.compound" => "C00139", "biocyc" => "META:Oxidized-ferredoxins", "seed.compound" => "cpd11621"))
GSM.addmet!(model, met)
met = GSM.createmetabolite("fdxrd_h", "Reduced ferredoxin hydrogenosome", "h", 0, "FeS",  Dict("kegg.compound" => "C00138", "biocyc" => "META:Reduced-ferredoxins", "seed.compound" => "cpd11620"))
GSM.addmet!(model, met)
met = GSM.createmetabolite("rq_ox_h", "Oxidized electron carrier hydrogenosome", "h", 0, "X")
GSM.addmet!(model, met)
met = GSM.createmetabolite("rq_red_h", "Reduced electron carrier hydrogenosome", "h", 0, "XH")
GSM.addmet!(model, met)

# Transporters
md = Dict("h2_c" => -1, "h2_h" => 1)
rxn = GSM.createreaction("H2ht", "H2 hydrogenosome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)
md = Dict("co2_c" => -1, "co2_h" => 1)
rxn = GSM.createreaction("CO2ht", "CO2 hydrogenosome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)
md = Dict("pi_c" => -1, "pi_h" => 1, "h_c" => -1, "h_h" => 1)
rxn = GSM.createreaction("PIht", "Phosphate hydrogenosome symporter", md)
GSM.addrxn!(model, rxn)
md = Dict("h2o_c" => -1, "h2o_h" => 1)
rxn = GSM.createreaction("H2Oht", "H2O hydrogenosome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)
md = Dict("adp_c" => -1, "adp_h" => 1, "atp_c" => 1, "atp_h" => -1)
rxn = GSM.createreaction("ADPATPht", "ADP/ATP translocase", md)
GSM.addrxn!(model, rxn)
md = Dict("succ_c" => -1, "pi_h" => -1, "pi_c" => 1, "succ_h" => 1)
rxn = GSM.createreaction("SUCCht1", "Succinate hydrogenosome transporter", md)
GSM.addrxn!(model, rxn)
md = Dict("succ_c" => -1, "h_c" => -1, "succ_h" => 1, "h_h" => 1)
rxn = GSM.createreaction("SUCCht2", "Succinate hydrogenosome transporter", md)
GSM.addrxn!(model, rxn)
md = Dict("mal__L_c" => -1, "pi_c" =>1, "pi_h" => -1, "mal__L_h" => 1)
rxn = GSM.createreaction("MALht", "Malate hydrogenosome antiporter", md)
GSM.addrxn!(model, rxn)

md = Dict("for_c" => -1, "for_h" => 1)
rxn = GSM.createreaction("FORht", "Formate hydrogenosome transporter", md)
GSM.addrxn!(model, rxn)
md = Dict("for_c" => -1, "h_c" => -1, "for_h" => 1, "h_h" => 1)
rxn = GSM.createreaction("FORht2", "Formate hydrogenosome symport", md)
GSM.addrxn!(model, rxn)

md = Dict("ac_c" => -1, "ac_h" => 1)
rxn = GSM.createreaction("ACht", "Acetate hydrogenosome transporter", md)
GSM.addrxn!(model, rxn)
md = Dict("ac_c" => -1, "h_c" => -1, "ac_h" => 1, "h_h" =>1)
rxn = GSM.createreaction("ACht2", "Acetate hydrogenosome symport", md)
GSM.addrxn!(model, rxn)

md = Dict("pyr_c" => -1, "h_c" => -1, "pyr_h" => 1, "h_h" =>1)
rxn = GSM.createreaction("PYRht", "Pyruvate hydrogenosome symport", md)
GSM.addrxn!(model, rxn)

# # Reactions
md = Dict("pyr_h" => -1, "coa_h" => -1, "accoa_h" => 1, "for_h" => 1) # PFL
rxn = GSM.createreaction("PFLh", "Pyruvate formate lyase", md)
GSM.addrxn!(model, rxn)

md = Dict("pyr_h" => -1, "coa_h" => -1, "fdxox_h" => -2, "accoa_h" => 1, "fdxrd_h" => 2, "co2_h" => 1, "h_h" => 1)
rxn = GSM.createreaction("PFOh", "Pyruvate ferredoxin oxidoreductase", md)
GSM.addrxn!(model, rxn)

#################################################################
# Hydrogenase: Option 1
md = Dict("h_h" => -2, "fdxrd_h" => -2, "fdxox_h" => 2, "h2_h" => 1)
rxn = GSM.createreaction("HYDhfe", "Ferredoxin Hydrogenase", md)
GSM.addrxn!(model, rxn)

md = Dict("nadh_h" => -1, "h_h" => -5, "rq_ox_h" => -2, "rq_red_h" => 2, "h_c" => 4, "nad_h" => 1) # rhodoquinone
rxn = GSM.createreaction("CMPL1h", "Complex 1 hydrogenosome", md)
GSM.addrxn!(model, rxn)

md = Dict("rq_ox_h" => 2, "rq_red_h" => -2, "fum_h" => -1, "succ_h" => 1) # H's come from rq (FUM:: C4H2O4 vs SUCC:: C4H4O4)
rxn = GSM.createreaction("CMPL2h", "Complex 2 hydrogenosome", md)
GSM.addrxn!(model, rxn)

md = Dict("fum_h" => -1, "h2o_h" => -1, "mal__L_h" => 1)
rxn = GSM.createreaction("FUMh", "Fumarase hydrogenosome", md)
GSM.addrxn!(model, rxn)

#################################################################
# Hydrogenase: Option 2
md = Dict("nadh_h"=>-1, "h_h"=>-3, "h2_h"=>2, "nad_h"=>1, "fdxrd_h" => -2, "fdxox_h" => 2) # NB only one hydrogenase
rxn = GSM.createreaction("HYDhbi", "Bifurcating Hydrogenase", md)
GSM.addrxn!(model, rxn)
#################################################################
# Option 3 - no complex 1 and 2
md = Dict("h_h" => -2, "fdxrd_h" => -2, "fdxox_h" => 2, "h2_h" => 1)
rxn = GSM.createreaction("HYDhfe", "Ferredoxin Hydrogenase", md)
GSM.addrxn!(model, rxn)
#################################################################

#################################################################
# ATP synthase?
md = Dict("atp_h" => 1, "h2o_h" => 1, "h_h" => 3, "adp_h" => -1, "h_c" => -4, "pi_h" => -1)
rxn = GSM.createreaction("ATPShydr", "ATP synthase hydrogenosome", md)
GSM.addrxn!(model, rxn)
#################################################################


md = Dict("adp_h" => -1, "pi_h" => -1, "succoa_h" => -1, "atp_h" => 1, "succ_h" => 1, "coa_h" => 1)
rxn = GSM.createreaction("SUCCOAh", "Succinyl-CoA synthase", md)
GSM.addrxn!(model, rxn)
#
md = Dict("succoa_h" => 1, "ac_h" => 1,  "succ_h" => -1, "accoa_h" => -1)
rxn = GSM.createreaction("SUCOAACTh", "Succinyl-CoA:acetate CoA transferase", md)
GSM.addrxn!(model, rxn)

md = Dict("mal__L_h" => -1, "nadh_h" => 1, "nad_h" => -1, "pyr_h" => 1, "co2_h" => 1)
rxn = GSM.createreaction("MEhy", "Malic enzyme (NAD) hydrogenosome", md)
GSM.addrxn!(model, rxn)

md = Dict("mal__L_h" => -1, "nadph_h" => 1, "nadp_h" => -1, "pyr_h" => 1, "co2_h" => 1)
rxn = GSM.createreaction("MEhx", "Malic enzyme (NADP) hydrogenosome", md)
GSM.addrxn!(model, rxn)


##################################################################
# # Shuttle
# md = Dict("mal__L_h" => -1, "oaa_h" => 1, "nad_h" => -1, "h_h" => 1, "nadh_h" => 1)
# rxn = GSM.createreaction("MDHh", "Malate dehydrogenase hydrogenosome", md)
# GSM.addrxn!(model, rxn)
#
# md = Dict("akg_h" => -1, "asp__L_h" => -1, "oaa_h" => 1, "glu__L_h" => 1)
# rxn = GSM.createreaction("AMTh", "Aspartate aminotransferase hydrogenosome", md)
# GSM.addrxn!(model, rxn)
#
# # Specific transporters
# md = Dict("mal__L_c" => -1, "mal__L_h" => 1, "akg_c" => 1, "akg_h" => -1)
# rxn = GSM.createreaction("MAKht", "Malate, alpha-ketoglutarate antiporter", md)
# GSM.addrxn!(model, rxn)
# md = Dict("asp__L_c" => -1, "asp__L_h" => 1, "glu__L_c" => 1, "glu__L_h" => -1)
# rxn = GSM.createreaction("AGht", "Aspartate, glutamate antiporter", md)
# GSM.addrxn!(model, rxn)

# Combined shuttle
md = Dict("nad_h" =>-1, "h_c" => -1, "nadh_c" => -1, "nadh_h" => 1, "h_h" => 1, "nad_c" => 1)
rxn = GSM.createreaction("MASh", "Mal/Asp shuttle", md)
GSM.addrxn!(model, rxn)
##################################################################

# basic hydrogenase that only goes in one direction reaction
md = Dict("h_h" => -1, "nadh_h" => -1, "nad_h" => 1, "h2_h" => 1)
rxn = GSM.createreaction("HYDh2", "Hydrogenase hydrogenosome", md)
GSM.addrxn!(model, rxn)

md = Dict("h_h" => -1, "nadph_h" => -1, "nadp_h" => 1, "h2_h" => 1)
rxn = GSM.createreaction("HYDh3", "Hydrogenase hydrogenosome", md)
GSM.addrxn!(model, rxn)
