require(Hotgenes)



# UpdateLevelsbyList ------------------------------------------------------

# change these levels
named_refs <- list(
  Type = c("Mississippi"),
  Treatment = c("chilled")
)

New_CO2 <- CO2 %>%
  UpdateLevelsbyList(named_refs)

m_2  <- model.matrix(~ Treatment * Type, data = New_CO2)

head(m_2 )


# valid_factors -----------------------------------------------------------

df_stws <- dplyr::starwars %>%
  dplyr::mutate_if(is.character, as.factor) %>%
  valid_factors()

df_stws

# these wil be dropped
df_stws$Invalid <- factor("One")
df_stws$Inv_NA <- NA

df_stws %>%
  valid_factors()
