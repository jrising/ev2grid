import pandas as pd
import os

data_path = os.path.join(os.path.dirname(__file__), "../data/predprice.csv")
out_path  = os.path.join(os.path.dirname(__file__), "../data/avg_regprice_byhour.csv")

df = pd.read_csv(data_path, parse_dates=["datetime"])
df["hour"] = df["datetime"].dt.hour

avg = df.groupby("hour")["predpe"].mean().reset_index()
avg.columns = ["hour", "avg_predpe"]

avg.to_csv(out_path, index=False)
print(f"Saved average regulation price by hour to {out_path}")
print(avg.to_string(index=False))
