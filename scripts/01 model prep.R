library(RMark)
library(R2ucare)
library(dplyr)


# overall CJS test

rekn <- read_inp('raw_data/SC eh cleaned.inp')
rekn_gof <- overall_CJS(X = rekn$encounter_histories, freq = rekn$sample_size)
rekn_gof

# Compare focal POPAN models

rekn_df <- as.data.frame(rekn$encounter_histories)
rekn_df <- rekn_df %>% 
  transmute(ch = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13))

run.rekn=function()
{
  #
  # Process data
  #
  rekn.processed=process.data(rekn_df, begin.time = 1, model = "POPAN")
  #
  # Create default design data
  #
  rekn.ddl=make.design.data(rekn.processed)
  #
  #
  #  Define range of models for Phi
  #
  Phi.dot=list(formula=~1)
  Phi.time=list(formula=~time)
  #
  #  Define range of models for p
  #
  p.dot=list(formula=~1)
  p.time=list(formula=~time)
  #
  #  Define range of models for pent
  #
  pent.time=list(formula=~time)
  #
  # Run all pairings of models
  #
  rekn.model.list=create.model.list("POPAN")
  rekn.results=mark.wrapper(rekn.model.list,
                              data=rekn.processed,ddl=rekn.ddl)
  #
  # Return model table and list of models
  #
  return(rekn.results)
}
rekn.results=run.rekn()

rekn.results


