# Read processed data from Ni et al
raw_TNT_table <- read.csv(
  './../raw_craniometric_data_straight_from_tnt.csv')#  #This data is processes raw data from Ni et al. 2021

raw_TNT_table$specimen <- sub("^[^_]*_", "", raw_TNT_table$specimen)

replace_sepcimen_name <- function(df, from, to) {
  df$specimen[df$specimen == from] <- to
  return(df)
}

raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._xuchang", to = "xuchang")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "dental", to = "Antecessor_dental")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "turkana", to = "knm_wt_15000")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._narmada", to = "narmada")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._eliye_springs", to = "eliye_springs")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._ndutu", to = "ndutu")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._rabat", to = "Homo_sp._rabat")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._stw53", to = "Homo_sp._stw53")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._tabun_c2", to = "tabun_c2")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "type", to = "Neanderthal_type")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._xiahe", to = "xiahe")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "altaiensis", to = "Homo_altaiensis")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._dali", to = "dali")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._hualongdong", to = "hualongdong")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "longi", to = "Homo_longi")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._jinniushan", to = "jinniushan")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._maba", to = "maba")
raw_TNT_table <- replace_sepcimen_name(raw_TNT_table, from = "sp._xuchang", to = "xuchang")


replace_column_name <- function(df, from, to) {
  colnames(df)[colnames(df) == from] <- to
  return(df)
}

raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M1..GOL..g.op..Maximum.cranial.lengt", to = "GOL")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M1d..NOL..n.op..Nasio.occipital.length", to = "NOL")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M5..BNL..Basion.nasion.length", to = "BNL")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M7..FOL..Foramen.magnum.length..ba.o.", to = "FOL")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M8..XCB..Maximum.cranial.breadth", to = "XCB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M10..XFB..Maximum.frontal.breadth..Frontal", to = "XFB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M10b..STB..Bistephanic.breadth..st.st", to = "STB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M11..AUB..Biauricular.breadth", to = "AUB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M12..ASB..Biasterion.breadth..ast.ast...Temporal", to = "ASB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M14..WCB..Minimum.cranial.breadth", to = "WCB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M17..BBH..Basion.Bregma.height", to = "BBH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M13a..MDB..Mastoid.width..Temporal", to = "MDB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M19a..MDH..Mastoid.height..Temporal", to = "MDH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M29..FRC..n.b..Frontal.sagital.chord", to = "FRC")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "FRF..Nasion.subtense.fraction..Frontal", to = "FRF")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M30..PAC..Parietal.sagital.chord", to = "PAC")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M31..OCC..l.o..Occipital.sagital.chord", to = "OCC")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M32.5...FRA..Frontal.angle..b.m.n..degree", to = "FRA")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M33d..OCA..Occipital.angle..degree.", to = "OCA")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M33e..PAA..Parietal.angle..degree.", to = "PAA")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M40..BPL..Basion.prosthion.length..Cranial.vault", to = "BPL")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M43a..FMB..Bifrontal.breadth..Frontal", to = "FMB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M43b..NAS..Nasio.frontal.subtense", to = "NAS")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M44..EKB..Biorbital.breadth..ek.ek", to = "EKB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M45..ZYB..Bizygomatic.breadth..zy.zy", to = "ZYB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M45.1...JUB..Bijugal.breadth", to = "JUB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M46b..ZMB..Bimaxillary.breadth", to = "ZMB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M48..NPH..Upper.facial.height..Nasion.prosthion.height..n.pr", to = "NPH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M48d..WMH..Cheek.height", to = "WMH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M49a..DKB..Interorbital.breadth..d.d", to = "DKB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M51a..OBB..Orbital.breadth", to = "OBB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M52..OBH..Orbital.height", to = "OBH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M54..NBL..Nasal.breadth", to = "NLB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M55..NH..Nasal.height..n.ns", to = "NLH")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M61..MAB..Maxilloalveolar.breadth", to = "MAB")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M76a..SSA..Zygomaxillary.angle..degree.", to = "SSA")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M77a..NFA..Nasio.frontal.angle..fm.a.n.fm.a", to = "NFA")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M41c..XML..Maximum.malar.lengt", to = "XML")
raw_TNT_table <- replace_column_name(raw_TNT_table, from = "M41d..MLS..Malar.subtens", to = "MLS")




# Read processed data from Ni et al
empty_matrix <- read.csv(
  './../empty_matrix.csv') #This data is an empty df with  metdata and colnaems used in this study

# Step 1: Preserve metadata rows (first two rows of empty_matrix)
metadata_rows <- empty_matrix[1:2, ]

# Step 2: Replace columns 2:5 in raw_TNT_table with those from empty_matrix
# Extract relevant columns from empty_matrix
empty_cols_2_to_5 <- empty_matrix[, c("specimen", colnames(empty_matrix)[2:5])]

# Merge to bring in updated cols 2:5 into raw_TNT_table
updated_df <- merge(
  raw_TNT_table,
  empty_cols_2_to_5,
  by = "specimen",
  suffixes = c("", ".new")
)

# Replace columns 2:5 with values from empty_matrix
updated_df[, colnames(empty_matrix)[2:5]] <- updated_df[, colnames(empty_matrix)[2:5]]

# Remove extra .new columns
updated_df <- updated_df[, !grepl(".new", colnames(updated_df), fixed = TRUE)]

# Step 3: Add metadata rows back on top
final_df <- rbind(metadata_rows, updated_df)

# Save
write.csv(final_df, file = "./../data_after_initial_processing.csv", row.names = FALSE, na = "")