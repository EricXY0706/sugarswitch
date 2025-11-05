import pandas as pd

proteins = ["binder2", "binder3", "binder135", "binder138"]
for p in proteins:
    df = pd.read_csv(f"/sdata2/WORK/Xuyi/{p}/{p}_single_points.csv")
    df.columns = ["BordaScore", "Site", "SS", "Mutation", "ConservationScore", "CouplingScore",
                  "SASA_i", "SASA_i+1", "SASA_i+2", "SASA_(i-1:i+1)", "SASA_(i:i+2)",
                  "ddG", "dTm", "ddG_NXS", "dTm_NXS", "ddG_NXT", "dTm_NXT", "MutScore", "MutScore_NXS", "MutScore_NXT", "Clash"]
    df.to_csv(f"/sdata2/WORK/Xuyi/{p}/{p}_single_points.csv", index=False)
