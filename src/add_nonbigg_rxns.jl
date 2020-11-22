# Add custom reactions to the model
met = Dict()
met["id"] = "xylan_e"
met["name"] = "Xylan"
met["charge"] = 0
met["formula"] = "C17H34O17"
anno = Dict()
anno["kegg.compound"] = "C00707; G00279; G00963"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "e"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Xylan exchange")
rxn["metabolites"] = Dict("xylan_e" => -1)
rxn["lower_bound"] = 0.0
rxn["upper_bound"] = 0.0
rxn["id"] = "EX_xylan_e"
rxn["annotation"] = Dict("sbo" => "SBO:0000627")
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict()
GSM.addrxn!(model, rxn)

rxn = Dict()
rxn["name"] = string("Xylanase")
rxn["metabolites"] = Dict("xylan_e" => -1, "glc__D_e" => 2, "xyl__D_e" => 1)
rxn["lower_bound"] = 0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_xylanase"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.2.1.18"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("Notes" => "Custom reaction")
GSM.addrxn!(model, rxn)

rxn = Dict()
rxn["name"] = string("Dextrin rxn")
rxn["metabolites"] = Dict("dextrin_e" => -1, "h2o_e" => -2, "glc__D_e" => 2)
rxn["lower_bound"] = 0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_dextrinase"
rxn["annotation"] = Dict()
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "TMP")
GSM.addrxn!(model, rxn)

rxn = Dict()
rxn["name"] = string("Glyceraldehyde-3-phosphate dehydrogenase (NADP+)")
rxn["metabolites"] = Dict("g3p_c" => -1, "nadp_c" => -1, "h2o_c" => -1, "nadph_c" => 1, "h_c" => 2, "3pg_c" => 1)
rxn["lower_bound"] = 0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_gapd"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["1.2.1.9"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC1.2.1.9")
GSM.addrxn!(model, rxn)
#################################################################
### Lysine
met = Dict()
met["id"] = "hicit_c"
met["name"] = "hicit_c"
met["charge"] = -3
met["formula"] = "C7H7O7"
met["notes"] = Dict()
met["annotation"] = Dict("inchikey" => "OEJZZCGRGVFWHK-WVZVXSGGSA-K", "seed.compound" => "cpd03372", "kegg.compound" => "C05662")
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Homoaconitate hydratase")
rxn["metabolites"] = Dict("hcit_c" => -1, "hicit_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_HAH"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["4.2.1.36"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC4.2.1.36")
GSM.addrxn!(model, rxn)

rxn = Dict()
rxn["name"] = string("Homoisocitrate dehydrogenase")
rxn["metabolites"] = Dict("hicit_c" => -1, "nad_c" => -1, "co2_c" => 1, "nadh_c" => 1, "2oxoadp_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_HICITD"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["1.1.1.87"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC1.1.1.87")
GSM.addrxn!(model, rxn)
# #################################################################

# #################################################################
# ### Cysteine and Methionine
rxn = Dict()
rxn["name"] = string("Methionine S-methyltransferase")
rxn["metabolites"] = Dict("amet_c" => -1, "met__L_c" => -1, "ahcys_c" => 1, "mmet_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_MetSmethtrans"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["2.1.1.12"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC2.1.1.12")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "s2ab_c"
met["name"] = "(S)-2-Aminobutanoate"
met["charge"] = 0
met["formula"] = "C4H9NO2"
anno = Dict()
anno["kegg.compound"] = "C02357"
anno["biocyc"] = "CPD-256"
anno["inchi"] = "QWCKQJZIFLGMSD-VKHMYHEASA-N"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

# (S)-2-Aminobutanoate + 2-Oxoglutarate <=> 2-Oxobutanoate + L-Glutamate
rxn = Dict()
rxn["name"] = string("2-aminobutanoate:2-oxoglutarate aminotransferase")
rxn["metabolites"] = Dict("s2ab_c"=>-1, "akg_c"=>-1, "2obut_c"=>1, "glu__L_c"=>1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_2abut2akgt"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["2.6.1.42"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC2.6.1.42")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "glglutl2abut_c"
met["charge"] = 1
met["name"] = "gamma-L-Glutamyl-L-2-aminobutyrate"
met["notes"] = Dict()
met["formula"] = "C4H8NO"
met["annotation"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

# ATP + L-Glutamate + (S)-2-Aminobutanoate <=> ADP + Orthophosphate + gamma-L-Glutamyl-L-2-aminobutyrate
rxn = Dict()
rxn["name"] = string("L-glutamate:2-aminobutanoate gamma-ligase (ADP-forming)")
rxn["metabolites"] = Dict("atp_c" =>-1, "s2ab_c"=>-1, "adp_c" => 1, "pi_c" => 1, "glglutl2abut_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_glutabutl"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["6.3.2.2"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC6.3.2.2")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "opth_c"
met["name"] = "Ophthalmate"
met["charge"] = 2
met["formula"] = "C4H7N"
anno = Dict()
anno["biocyc"] = "CPD-20340"
anno["inchi"] = "JCMUOFQHZLPHQP-BQBZGAKWSA-M"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

# ATP + gamma-L-Glutamyl-L-2-aminobutyrate + Glycine <=> ADP + Orthophosphate + Ophthalmate
rxn = Dict()
rxn["name"] = string("gamma-L-glutamyl-L-2-aminobutyrate:glycine ligase (ADP-forming)")
rxn["metabolites"] = Dict("glglutl2abut_c" => -1, "atp_c" =>-1, "adp_c" => 1, "pi_c" => 1, "opth_c"=>1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_glglutglyl"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["6.3.2.3"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC6.3.2.3")
GSM.addrxn!(model, rxn)
# #################################################################
#
# #################################################################
# ### Starch and sucrose
met = Dict()
met["id"] = "cellulose_e"
met["name"] = "Cellulose"
met["charge"] = 0
met["formula"] = "C24H44O22"
anno = Dict()
anno["Name"] = "Name: (1,4-beta-D-Glucosyl)n; (1,4-beta-D-Glucosyl)n+1; (1,4-beta-D-Glucosyl)n-1; 1,4-beta-D-Glucan; 1,4-beta-D-glucan; 1,4-beta-Glucan; Cellulose; Microcrystalline cellulose; beta-D-glucan; beta-glucan; glucan"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "e"
GSM.addmet!(model, met)

met = Dict()
met["id"] = "hemicellulose_e"
met["name"] = "Hemicellulose"
met["charge"] = 0
met["formula"] = "C10H20O10"
anno = Dict()
anno["Name"] = "Hemicellulose"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "e"
GSM.addmet!(model, met)


rxn = Dict()
rxn["name"] = string("Cellulose exchange")
rxn["metabolites"] = Dict("cellulose_e" => -1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "EX_cellulose_e"
rxn["annotation"] = Dict("sbo" => "SBO:0000627")
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict()
GSM.addrxn!(model, rxn)

rxn = Dict()
rxn["name"] = string("Hemicellulose exchange")
rxn["metabolites"] = Dict("hemicellulose_e" => -1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "EX_hemicellulose_e"
rxn["annotation"] = Dict("sbo" => "SBO:0000627")
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict()
GSM.addrxn!(model, rxn)


# Cellulose_e => cellobiose_e
rxn = Dict()
rxn["name"] = string("Cellulose to cellobiose lumped reaction EXTERNAL")
rxn["metabolites"] = Dict("cellulose_e" => -1, "cellb_e" => 2)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_cellulase"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.2.1.4"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.2.1.4")
GSM.addrxn!(model, rxn)

# Hemicellulose_e => xylose_e
rxn = Dict()
rxn["name"] = string("Hemicellulose to xylose lumped reaction EXTERNAL")
rxn["metabolites"] = Dict("hemicellulose_e" => -1, "xyl__D_e" => 2)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_hemicellulase"
rxn["annotation"] = Dict{String, Any}("ec-code"=> ["3.2.1.3","3.2.1.6","3.2.1.8","3.2.1.10","3.2.1.11","3.2.1.20","3.2.1.22","3.2.1.23","3.2.1.24","3.2.1.25","3.2.1.26","3.2.1.32","3.2.1.33","3.2.1.37","3.2.1.39","3.2.1.40","3.2.1.43","3.2.1.18"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "Assorted")
GSM.addrxn!(model, rxn)
# #################################################################
#
# #################################################################
# # Pyruvate metabolism
met = unimets["lac__L_c"]
GSM.addmet!(model, met)

# # lac__D_c + 2.0 ficytb5_c ⇌ pyr_c + 2.0 focytb5_c
rxn = Dict()
rxn["name"] = string("D-Lactate dehydrogenase using heme")
rxn["metabolites"] = Dict("lac__D_c" => -1,  "ficytb5_c" => -2, "pyr_c" => 1,  "focytb5_c" => 2, "h_c" => 2)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_lacd1"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["1.1.2.4"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC1.1.2.4")
GSM.addrxn!(model, rxn)
#
# # lac__L_c + 2.0 ficytb5_c ⇌ pyr_c + 2.0 focytb5_c
rxn = Dict()
rxn["name"] = string("L-lactate dehydrogenase using heme")
rxn["metabolites"] = Dict("lac__L_c" => -1,  "ficytb5_c" => -2, "pyr_c" => 1,  "focytb5_c" => 2, "h_c" => 2)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_lacd2"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["1.1.2.4"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC1.1.2.4")
GSM.addrxn!(model, rxn)
# #################################################################
#
# #################################################################
# # Sphingolipid
met = Dict()
met["id"] = "acoa_c"
met["name"] = "Acyl-CoA"
met["notes"] = Dict()
met["charge"] = -4
met["formula"] = "C22H31N7O17P3SR"
anno = Dict("kegg.compound" => "C00040", "seed.compound" => "cpd11611", "biocyc" => "META:ACYL-COA", "chebi" => "CHEBI:13727, CHEBI:13802, CHEBI:17984, CHEBI:22223, CHEBI:24025, CHEBI:2455, CHEBI:37554, CHEBI:4987, CHEBI:58342, CHEBI:77636")
met["compartment"] = "c"
GSM.addmet!(model, met)

met = Dict()
met["id"] = "cer_18_p_c"
met["name"] = "Ceramide 1 Phosphate"
met["charge"] = -2
met["formula"] = "C36H70NO6P"
anno = Dict()
anno["bigg.metabolite"] = "crmp_hs"
anno["kegg.compound"] = "C02960"
anno["biocyc"] = "CPD-502"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("3-sn-phosphatidate phosphohydrolase")
rxn["metabolites"] = Dict("cer_18_c" =>-1, "pi_c" =>-1, "h2o_c" => 1, "cer_18_p_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R06522"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.1.3.4"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.1.3.4")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "sphmye_c"
met["name"] = "Sphingomyelin"
met["charge"] = -1
met["formula"] = "C36H71NO6PR"
anno = Dict()
anno["biocyc"] = "N-acyl-sphingosylphosphorylcholine"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)


met = Dict()
met["id"] = "ppchol_c"
met["name"] = "Phosphocholine"
met["charge"] = -1
met["formula"] = "H2O4PR"
anno = Dict()
anno["inchi"] = "YHHSONZFOIEMCP-UHFFFAOYSA-M"
anno["kegg.compound"] = "C00588"
anno["biocyc"] = "PHOSPHORYL-CHOLINE"
anno["bigg.metabolite"] = "cholp"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)


rxn = Dict()
rxn["name"] = string("Sphingomyelin phosphodiesterase")
rxn["metabolites"] = Dict("cer_18_c" =>-1, "ppchol_c" =>-1, "h2o_c" => 1, "sphmye_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_02541"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.1.3.4"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.1.3.4")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "spppchol_c"
met["name"] = "Sphingosyl-phosphocholine"
met["charge"] = -1
met["formula"] = "C35H72NO5P"
anno = Dict()
anno["inchi"] = "JLVSPVFPBBFMBE-HXSWCURESA-O"
anno["biocyc"] = "CPD-481"
anno["kegg.compound"] = "C03640"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Acyl-CoA:sphingosine N-acyltransferase")
rxn["metabolites"] = Dict("spppchol_c" =>1, "acoa_c" =>1, "coa_c" => -1, "sphmye_c" => -1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R02543"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["2.3.1.24"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC2.3.1.24")
GSM.addrxn!(model, rxn)


met = Dict()
met["id"] = "gluccera_c"
met["name"] = "D-glucosyl-N-acylsphingosine"
met["charge"] = 0
met["formula"] = "C42H81NO8"
anno = Dict()
anno["kegg.compound"] = "C01190; G10238"
anno["biocyc"] = "Glucosyl-acyl-sphingosines"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Glucosylceramidase")
rxn["metabolites"] = Dict("cer_18_c" =>-1, "glc__D_c" =>-1, "gluccera_c" => 1, "h2o_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R01498"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.2.1.45"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.2.1.45")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "lacsylcera_c"
met["name"] = "Lactosylceramide"
met["charge"] = 0
met["formula"] = "C48H91NO13"
anno = Dict()
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Beta-galactosidase")
rxn["metabolites"] = Dict("gluccera_c" =>-1, "gal_c" =>-1, "lacsylcera_c" => 1, "h2o_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R01100"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.2.1.23"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.2.1.23")
GSM.addrxn!(model, rxn)


met = Dict()
met["id"] = "galsylcera_c"
met["name"] = "Galactosylceramide"
met["charge"] = 1
met["formula"] = "C42H82NO8"
anno = Dict()
anno["inchi"] = "DTFAJAKTSMLKAT-JDCCYXBGSA-P"
anno["kegg.compound"] = "C02686; G11121"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Beta-galactosidase")
rxn["metabolites"] = Dict("cer_18_c" =>-1, "udpgal_c" =>-1, "galsylcera_c" => 1, "udp_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R01500"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["2.4.1.47"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC2.4.1.47")
GSM.addrxn!(model, rxn)

met = Dict()
met["id"] = "digalsylcera_c"
met["name"] = "Digalactosylceramide"
met["charge"] = 1
met["formula"] = "C48H92NO13"
anno = Dict()
anno["kegg.compound"] = "C06126"
anno["BIGG"] = "digalside_hs"
anno["biocyc"] = "Digalactosylceramides"
met["annotation"] = anno
met["notes"] = Dict()
met["compartment"] = "c"
GSM.addmet!(model, met)

rxn = Dict()
rxn["name"] = string("Alpha-galactosidase")
rxn["metabolites"] = Dict("galsylcera_c" =>-1, "gal_c" =>-1, "digalsylcera_c" => 1, "h2o_c" => 1)
rxn["lower_bound"] = -1000.0
rxn["upper_bound"] = 1000.0
rxn["id"] = "r_R04019"
rxn["annotation"] = Dict{String, Any}("ec-code" => ["3.2.1.22"])
rxn["gene_reaction_rule"] = ""
rxn["notes"] = Dict("EC" => "EC3.2.1.22")
GSM.addrxn!(model, rxn)


#################################################################
# Vitamin requirements
GSM.addsupply!(model, "4mhetz_c") # need to supply yhis
GSM.addsink!(model, "dad_5_c")
# GSM.addsupply!(model, "malACP_c")
GSM.addsink!(model, "amob_c")
GSM.addsupply!(model, "4ahmmp_c") # need to supply this
#################################################################
