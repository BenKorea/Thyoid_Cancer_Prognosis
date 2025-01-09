################################################################################
## my_functions_for_common_purpose
################################################################################
my_clean_string_edges <- function(str) {
  
  # 정규식 설명:
  # ^\\s*          : 문자열의 시작과 0개 이상의 공백 문자
  # (?:#1)?        : '#1' 패턴이 0번 또는 1번 나타나는 부분 (캡처를 하지 않는 그룹)
  # [.,]?          : '.' 또는 ','가 0번 또는 1번 나타나는 부분
  # \\s*           : 다시 0개 이상의 공백 문자
  str <- sub("^\\s*(?:#1)?[.,]?\\s*", "", str)
  str <- gsub("\\s+$", "", str)
  
  return(str)
}


################################################################################
## my_functions_for_hash
################################################################################
my_deidentify_PatientID <- function(df, col_name) {
  # 열을 정렬하고 고유한 값들을 추출합니다.
  unique_values <- unique(df[[col_name]])
  sorted_values <- sort(unique_values)
  
  # 고유한 값들에 대해 순서를 부여합니다.
  value_to_index <- match(sorted_values, unique_values)
  
  # 해당 열을 고유한 값에 대응하는 숫자로 대체합니다.
  df[[col_name]] <- as.character(value_to_index[match(df[[col_name]], unique_values)])
  
  return(df)
}


################################################################################
## my_functions_for_birthday
################################################################################
my_parsing_birthday_data<-function(dt) {
 
  birthday_data <- dt[분류명 == "Birthday_Registration"]
  birthday_data[, birthday := as.Date(특기사항, format="%Y-%m-%d")]
  if (input_error_checking_mode == "Y") {
    birthday_input_error <<- birthday_data[is.na(birthday), ]
  }
  birthday_data <- birthday_data[!is.na(birthday), ]
  birthday_data[, 다운로드기준나이 := as.integer(difftime(downloaded_date, birthday, units = "weeks") / 52.25)]
  birthday_data[, validation := 나이-다운로드기준나이]
  if (input_error_checking_mode == "Y") {
    birthday_data_validation_error <<- birthday_data[validation > 1 | validation < -1]
  }
  birthday_data <- birthday_data[!(validation > 1 | validation < -1)]
  birthday_data <- birthday_data[, c('등록번호', '성별', 'birthday')]
  return(birthday_data)
  
}


################################################################################
## my_functions_for_risk
################################################################################
my_splilt_lines_for_risk_data <- function(dt) {

  dt[, line_count := str_count(특기사항, "\n") + 1]
  if (input_error_checking_mode == "Y") {
    risk_data_line_count_errors <<- dt[line_count >= 4 | line_count == 1]
  }
  dt <- dt[line_count == 3 | line_count == 2]
  dt[line_count == 3, c("cN_line", "op_line", "stage_line") := tstrsplit(특기사항, split = "\n", fixed = TRUE, fill = NA)]
  dt[line_count == 2, c("op_line", "stage_line") := tstrsplit(특기사항, split = "\n", fixed = TRUE, fill = NA)]
  
  dt[, cN_line := my_clean_string_edges(cN_line)]
  dt[, op_line := my_clean_string_edges(op_line)]
  dt[, stage_line := my_clean_string_edges(stage_line)]
  
  if (input_error_checking_mode == "Y") {
    risk_data_op_line_error <<- dt[!grepl("^[0-9]", dt$op_line), ]
  }
  dt <- dt[grepl("^[0-9]", dt$op_line), ]
  
  return(dt)
}

my_parsing_cN_line <- function(dt) {
  
  dt$cN_date <- tstrsplit(dt$cN_line, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$cN_modality <- tstrsplit(dt$cN_line, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$cN <- tstrsplit(dt$cN_line, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  # "PCNA (+)"가 있으면 dt$PCNA를 "Positiv"로 파생  
  dt$PCNA <- ifelse(grepl("PCNA \\(\\+\\)", dt$cN_line), "Positive", NA)
  # &가 있으면 이후를 PETCT_part로 파생
  dt$PETCT_part <- ifelse(grepl("&", dt$cN_line), sub(".*&\\s*", "", dt$cN_line), NA)
  dt$PETCT_date <- tstrsplit(dt$PETCT_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]  
  dt$Thyroid_SUVmax <- tstrsplit(dt$PETCT_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  dt$Thyroid_SUVmax <- sub("Thyroid=", "", dt$Thyroid_SUVmax)
  dt$Thyroid_SUVmax <- sub(",", "", dt$Thyroid_SUVmax)
  dt$PETCT_Metastasis <- tstrsplit(dt$PETCT_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[4]]
  dt$PETCT_Metastasis_SUVmax <- tstrsplit(dt$PETCT_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[5]]
  dt$PETCT_Metastasis<-paste0(dt$PETCT_Metastasis,".",dt$PETCT_Metastasis_SUVmax)
  dt$PETCT_Metastasis_SUVmax <- ifelse(grepl("=", dt$PETCT_Metastasis), sub(".*=\\s*", "", dt$PETCT_Metastasis), NA)
  dt$PETCT_Metastasis <- ifelse(grepl("=", dt$PETCT_Metastasis), sub("\\s*=.*", "", dt$PETCT_Metastasis), NA)
  
  dt[, cN_date := as.Date(cN_date, format="%Y-%m-%d", try = TRUE)]
  if (input_error_checking_mode == "Y") {
    risk_data_cN_date_error <<- dt[is.na(cN_date) & line_count==3]
  }
  dt <- dt[!(is.na(cN_date) & line_count==3)]
  
  if (input_error_checking_mode == "Y") {
    risk_data_cN_error <<- dt[!cN %in% c("cN0", "cN1a", "cN1b", "cN1", NA,"pN0", "pN1a", "pN1b", "pN1")]
  }
  dt <- dt[cN %in% c("cN0", "cN1a", "cN1b", "cN1", NA,"pN0", "pN1a", "pN1b", "pN1")]
  
  return(dt)
}

my_parsing_op_line <- function(dt) {
 
  # &를 기준으로 op_part와 completion 열을 분리합니다.
  dt[, c("op_part", "completion_part") := tstrsplit(op_line, "&", fixed = TRUE)]
  dt[, op_part := my_clean_string_edges(op_part)]
  dt[, completion_part := my_clean_string_edges(completion_part)]

  dt[, c("op", "op_hosp") := tstrsplit(op_part, "@", fixed = TRUE)]
  dt[, op_part := my_clean_string_edges(op_part)]
  dt[, op_hosp := my_clean_string_edges(op_hosp)]
  
  dt$op_date <- tstrsplit(dt$op_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$op_name <- tstrsplit(dt$op_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$operator <- tstrsplit(dt$op_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  dt[, op_date := my_clean_string_edges(op_date)]
  dt[, op_name := my_clean_string_edges(op_name)]
  dt[, operator := my_clean_string_edges(operator)]
  
  # op_date 열에서 "-??"가 존재하면 "-06"으로 대체
  dt$op_date <- ifelse(grepl("-\\?\\?", dt$op_date), gsub("-\\?\\?", "-06", dt$op_date), dt$op_date)
  dt[, op_date := as.Date(op_date, format="%Y-%m-%d", try = TRUE)]
  if (input_error_checking_mode == "Y") {
  risk_data_op_date_error <<- dt[is.na(op_date)]
  }
  dt <- dt[!is.na(op_date)]
  
  # 
  dt[, c("completion_part", "completion_hosp") := tstrsplit(completion_part, "@", fixed = TRUE)]
  dt[, completion_part := my_clean_string_edges(completion_part)]
  dt[, completion_hosp := my_clean_string_edges(completion_hosp)]
  
  dt$completion_date <- tstrsplit(dt$completion_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$completion_name <- tstrsplit(dt$completion_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$completion_operator <- tstrsplit(dt$completion_part, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  dt[, completion_date := my_clean_string_edges(completion_date)]
  dt[, completion_name := my_clean_string_edges(completion_name)]
  dt[, completion_operator := my_clean_string_edges(completion_operator)]
  
  dt[, completion_date := as.Date(completion_date, format="%Y-%m-%d", try = TRUE)]
  if (input_error_checking_mode == "Y") {
    risk_data_completion_date_error <<- dt[is.na(completion_date)&!is.na(completion_part)]
  }
  dt <- dt[!is.na(completion_date)|is.na(completion_part)]
  
  return(dt)
}

my_mutate_op_and_ND_type <- function(dt) {

  dt$Thyroidectomy_Type <- "Others"
  dt$Thyroidectomy_Type[grepl("sub", dt$op_name, ignore.case = TRUE)|grepl("near", dt$op_name, ignore.case = TRUE)] <- "subTT"
  dt$Thyroidectomy_Type[!grepl("sub", dt$op_name, ignore.case = TRUE)&!grepl("near", dt$op_name, ignore.case = TRUE)&grepl("TT", dt$op_name, ignore.case = TRUE)] <- "TT"
  dt$Thyroidectomy_Type[!is.na(dt$completion_name)] <- "Lobectomy_Completion"
  dt$Thyroidectomy_Type[is.na(dt$completion_name)&grepl("lobectomy", dt$op_name, ignore.case = TRUE)] <- "Lobectomy"
  dt$Thyroidectomy_Type[is.na(dt$completion_name)&grepl("completion", dt$op_name, ignore.case = TRUE)] <- "Completion"
  
  dt$Endoscopic_or_Robotic<-"Conventional_Surgery"  
  dt$Endoscopic_or_Robotic[grepl("E_", dt$op_name, ignore.case = TRUE)] <- "Endoscopic"
  dt$Endoscopic_or_Robotic[grepl("R_", dt$op_name, ignore.case = TRUE)] <- "Robotic"
  
  dt$ND_Type <- "Not Done"
  dt$ND_Type[!is.na(dt$completion_name)] <-"CND"
  dt$LND_extent<-tstrsplit(dt$op_name, "\\_", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$ND_Type[grepl("cND", dt$op_name, ignore.case = TRUE)] <- "CND"
  dt$ND_Type[grepl("ND", dt$LND_extent, ignore.case = TRUE)] <- "ND"
  dt$ND_Type[grepl("FND", dt$LND_extent, ignore.case = TRUE)] <- "FND"
  dt$ND_Type[grepl("SND", dt$LND_extent, ignore.case = TRUE)] <- "SND"
  dt$ND_Type[grepl("mRND", dt$LND_extent, ignore.case = TRUE)] <- "mRND"
  dt$ND_Type[grepl("picking", dt$LND_extent, ignore.case = TRUE)] <- "Picking"
  
  return(dt)
}

my_split_paracentesis_for_stage_line <- function(dt) {

  # ")"를 기준으로 stage_line 열 분리
  dt$histology_part <- tstrsplit(dt$stage_line, "\\) ", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$pT_part <- tstrsplit(dt$stage_line, "\\) ", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$pN_part <- tstrsplit(dt$stage_line, "\\) ", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  dt$M_part <- tstrsplit(dt$stage_line, "\\) ", fixed = FALSE, fill = NA, type.convert = TRUE)[[4]]
  
  return(dt)
  
}

my_parsing_histology <- function(dt) {
  
  # split into highest and subsequent stage histology
  dt$highest_stage_histology_part <- tstrsplit(dt$histology_part, "\\+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$subsequent_stage_histology_part <- tstrsplit(dt$histology_part, "\\+", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$subsequent_stage_histology_part<-ifelse(grepl("NIFTP",dt$subsequent_stage_histology_part, ignore.case = TRUE),NA,dt$subsequent_stage_histology_part)
  
  # split into histology and location and number
  dt$highest_stage_histology <- tstrsplit(dt$highest_stage_histology_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  # side and number
  dt$highest_stage_multiplicity <- ifelse(dt$Risk=="Unknown",NA,tstrsplit(dt$highest_stage_histology_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]])
  # 공백을 언더바로 대체
  dt$highest_stage_multiplicity <- ifelse(!grepl("_",dt$highest_stage_multiplicity, ignore.case = TRUE), sub(" ","_",dt$highest_stage_multiplicity),dt$highest_stage_multiplicity)
  
  dt$highest_1st <- tstrsplit(dt$highest_stage_multiplicity, "\\s", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$highest_1st_location <- tstrsplit(dt$highest_1st, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$highest_1st_number <- tstrsplit(dt$highest_1st, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$highest_1st_number <- gsub("x", "", dt$highest_1st_number)
  dt$highest_1st_number <- gsub(")", "", dt$highest_1st_number)
  
  
  dt$highest_2nd <- tstrsplit(dt$highest_stage_multiplicity, "\\s", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$highest_2nd_location <- tstrsplit(dt$highest_2nd, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$highest_2nd_number <- tstrsplit(dt$highest_2nd, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$highest_2nd_number <- gsub("x", "", dt$highest_2nd_number)
  dt$highest_2nd_number <- gsub(")", "", dt$highest_2nd_number)
  
    
  dt$highest_3rd <- tstrsplit(dt$highest_stage_multiplicity, "\\s", fixed = FALSE, fill = NA, type.convert = TRUE)[[3]]
  dt$highest_3rd_location <- tstrsplit(dt$highest_3rd, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$highest_3rd_number <- tstrsplit(dt$highest_3rd, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$highest_3rd_number <- gsub("x", "", dt$highest_3rd_number)
  dt$highest_3rd_number <- gsub(")", "", dt$highest_3rd_number)
  
  
  # split into histology and location and number
  dt$subsequent_stage_histology <- tstrsplit(dt$subsequent_stage_histology_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  # side and number
  dt$subsequent_stage_multiplicity <- ifelse(dt$Risk=="Unknown",NA,tstrsplit(dt$subsequent_stage_histology_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]])
  # 공백을 언더바로 대체
  dt$subsequent_stage_multiplicity <- ifelse(!grepl("_",dt$subsequent_stage_multiplicity, ignore.case = TRUE), sub(" ","_",dt$subsequent_stage_multiplicity),dt$subsequent_stage_multiplicity)
  
  dt$subsequent_1st <- tstrsplit(dt$subsequent_stage_multiplicity, "\\s", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$subsequent_1st_location <- tstrsplit(dt$subsequent_1st, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$subsequent_1st_number <- tstrsplit(dt$subsequent_1st, "_", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  dt$subsequent_1st_number <- gsub("x", "", dt$subsequent_1st_number)
 
  dt$Multiplicity<- ifelse(dt$highest_1st_number=="1"&is.na(dt$highest_2nd_number)&is.na(dt$highest_3rd_number),"Single","Multiple")
  dt$Laterality<- ifelse(grepl("Rt",dt$highest_stage_histology_part, ignore.case = TRUE)&grepl("Lt",dt$highest_stage_histology_part, ignore.case = TRUE),"Bilateral","Unilateral")
  dt$Laterality<- ifelse(grepl("Rt",dt$highest_stage_histology_part, ignore.case = TRUE)&grepl("Lt",dt$subsequent_stage_histology_part, ignore.case = TRUE),"Bilateral",dt$Laterality)
  dt$Laterality<- ifelse(grepl("Lt",dt$highest_stage_histology_part, ignore.case = TRUE)&grepl("Rt",dt$subsequent_stage_histology_part, ignore.case = TRUE),"Bilateral",dt$Laterality)
  dt$Laterality<- ifelse(grepl("both",dt$highest_stage_histology_part, ignore.case = TRUE)|grepl("both",dt$subsequent_stage_histology_part, ignore.case = TRUE),"Bilateral",dt$Laterality)
    
  return(dt)
  
}
  
my_parsing_pT <- function(dt) {  
  
  dt$pT <- tstrsplit(dt$pT_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$pT_detail <- tstrsplit(dt$pT_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  
  dt$pT_temp <- tstrsplit(dt$pT_detail, "\\s", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  
  # dt$pT_subtype, dt$pT_size 열 초기화
  dt$pT_subtype <- NA
  dt$pT_size <- NA
  
  # dt$pT_temp 열의 값이 알파벳으로 시작하면 dt$pT_subtype으로 파생, 숫자로 시작하면 dt$pT_size로 파생
  for (i in seq_along(dt$pT_temp)) {
    if (grepl("^[A-Za-z]", dt$pT_temp[i])) {  # 알파벳으로 시작하는지 확인
      dt$pT_subtype[i] <- dt$pT_temp[i]  # 알파벳으로 시작하면 dt$pT_subtype으로 파생
    } else {
      dt$pT_size[i] <- dt$pT_temp[i]  # 숫자로 시작하면 dt$pT_size로 파생
    }
  }
  
  dt$ETE <- "N"
  dt$ETE[grepl("ETE\\+", dt$pT_detail)] <- "Y"
  
  dt$gross_ETE <- "N"
  dt$gross_ETE[grepl("gross_ETE\\+", dt$pT_detail)] <- "Y"
  
  
  return(dt)
  
}

my_parsing_pN <- function(dt) {  
  
  dt$pN <- tstrsplit(dt$pN_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$pN_detail <- tstrsplit(dt$pN_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  
  dt$ENE <- "N"
  dt$ETE[grepl("ETE\\+", dt$pT_detail)] <- "Y"
  
  dt$ETE <- "N"
  dt$ENE[grepl("ENE+", dt$pN_detail)] <- "Y"
  return(dt)
  
}

my_parsing_M <- function(dt) {  
  
  dt$M <- tstrsplit(dt$M_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  dt$M_detail <- tstrsplit(dt$M_part, "\\(", fixed = FALSE, fill = NA, type.convert = TRUE)[[2]]
  
  dt$M_Bone <- "N"
  dt$M_Bone[grepl("Bone", dt$M_detail)] <- "Y"  
  
  dt$M_Lung <- "N"
  dt$M_Lung[grepl("Lung", dt$M_detail)] <- "Y"    
  # 
  # dt[, M := ifelse(is.null(M_part), "N", "Y")]
  return(dt)
  
}


my_parsing_risk_data<-function(dt) {
  
  risk_data <- dt[분류명 %in% c("Low", "Intermediate", "High", "Unknown")]
  setnames(risk_data, "분류명", "Risk")
  
  risk_data <- my_splilt_lines_for_risk_data(risk_data)
  risk_data <- my_parsing_cN_line(risk_data)
  risk_data <- my_parsing_op_line(risk_data)
  risk_data <- my_mutate_op_and_ND_type (risk_data)
  risk_data <- my_split_paracentesis_for_stage_line (risk_data)
  risk_data <- my_parsing_histology (risk_data)
  risk_data <- my_parsing_pT (risk_data)
  risk_data <- my_parsing_pN (risk_data)
  risk_data <- my_parsing_M (risk_data)
  
  return(risk_data)
}

################################################################################
## my_functions_for_response
################################################################################
my_parsing_response_line <- function(dt) {
  
  dt$response_date <- tstrsplit(dt$특기사항, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]

  return(dt)
}

my_mutate_recur_date <- function(dt) {
  
  dt$Recur<- ifelse(dt$Response=="Structural","Y","N")

  return(dt)
}

my_parsing_response_data<-function(dt) {
  
  response_data <- dt[분류명 %in% c("Excellent","Indeterminate","Biochemical","Structural")]
  setnames(response_data, "분류명", "Response")
  response_data[, 등록일 := as.Date(등록일, format="%Y-%m-%d")]
  response_data <- my_parsing_response_line(response_data)
  response_data <- my_mutate_recur_date(response_data)
  
  return(response_data)
  
}

################################################################################
## my_functions_for_followup
################################################################################
my_parsing_response_line <- function(dt) {
  
  dt$response_date <- tstrsplit(dt$특기사항, "\\s+", fixed = FALSE, fill = NA, type.convert = TRUE)[[1]]
  
  return(dt)
}

my_mutate_recur_date <- function(dt) {
  
  dt$Recur<- ifelse(dt$Response=="Structural","Y","N")
  
  return(dt)
}

my_parsing_followup_data<-function(dt) {
  
  followup_data <- dt[분류명 %in% c("Transfer","follow up loss","expire","요양병원","호스피스")]
  setnames(followup_data, "분류명", "FollowUp")
  followup_data[, 등록일 := as.Date(등록일, format="%Y-%m-%d")]

  
  return(followup_data)
  
}

  