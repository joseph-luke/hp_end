# %%
import pandas as pd
import re
from itertools import groupby
from pathlib import Path

# %%
data_folder = Path(
	"C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"
)

data = {}
with open(data_folder / "Processed/stylised_hp.txt", "r", encoding="utf-8") as file:
	for line in file:
		s = line.split(":", maxsplit=1)
		x = int(s[0])
		t = s[1].split("\n")[0].strip()
		data[x] = t
file.close()

file_to_open_1 = data_folder / "Processed/text_corrected_hp.csv"
df = pd.read_csv(file_to_open_1, low_memory=False)

# %%
abbreviated_categories = [
	"Cs",
	"Ca",
	"Sa",
	"Ly",
	"Le",
	"Bl",
	"Gl",
	"Ml",
	"Ms",
	"Te",
	"Se",
	"As",
	"PNST",
	"My",
	"∪",
]

# Remove records that do not mention malignancy in their main diagnostic portion
m_mentioning_records = []
no_m = []

for key, value in data.items():
	data_columns = []
	if re.search("MicroDescr", value):
		data_columns.append("MicroDescr")
	if re.search("Diagnosis", value):
		data_columns.append("Diagnosis")
	if re.search("Comments", value):
		data_columns.append("Comments")

	if not data_columns:
		no_m.append(key)

	# Obtain cellular data in the correct order
	data_cells = []
	if data_columns:
		full_row_data = re.findall(r"{(.*)}", value)[0]
		column_data_as_list = re.split(r",", full_row_data)

		for item in column_data_as_list:
			cell = re.findall(r": '(.*)'", item)
			for i in data_columns:
				if re.findall(i, item):
					data_cells.append(cell[0])

	# Now we have the columns and their respective information in 2 lists
	# Diagnostic portion == Diagnosis info => end (if diagnosis info exists),
	# else, Micro info => end
	# We filter out the records that do not mention a malignancy in their
	# diagnositic portion

	diagnostic = " ".join(data_cells).strip()

	n = 0  # Types of malignancy in row == 0
	for cat in abbreviated_categories:
		if re.search(cat, diagnostic):
			n += 1
	if n >= 1:
		m_mentioning_records.append(key)
	elif n < 1:
		no_m.append(key)


# %%
m_mentioning_records = list(sorted(set(m_mentioning_records)))

data1 = {}
for i in m_mentioning_records:
	data1[i] = data[i]
data = data1

# %%
# Re-number keys in data
df = df.loc[m_mentioning_records]

data1 = {}
new_keys = list(range(df.shape[0]))

for new_key, old_key in zip(new_keys, data.keys()):
	data1[new_key] = data[old_key]

data = data1

df.index = list(data.keys())

# %%
condensed_data = {}
condensed_key_subset = []

for key, value in data.items():
	data_columns = []
	if re.search("MicroDescr", value):
		data_columns.append("MicroDescr")
	if re.search("Diagnosis", value):
		data_columns.append("Diagnosis")
	if re.search("Comments", value):
		data_columns.append("Comments")

	# Obtain cellular data in the correct order
	data_cells = []
	if data_columns:
		full_row_data = re.findall(r"{(.*)}", value)[0]
		column_data_as_list = re.split(r",", full_row_data)

		for item in column_data_as_list:
			cell = re.findall(r": '(.*)'", item)
			for i in data_columns:
				if re.findall(i, item):
					data_cells.append(cell[0])

	# Select Micro: or Diagnosis: (Micro: if there is no Diagnosis info)
	if "Diagnosis" not in data_columns:  # Could have Micro or Comments or both
		# Join everything
		diagnostic = " ".join(data_cells).strip()

	elif "Diagnosis" in data_columns:
		if "MicroDescr" in data_columns:  # Get everything past Micro
			del data_columns[0]
			del data_cells[0]
			diagnostic = " ".join(data_cells).strip()

		elif "MicroDescr" not in data_columns:  # Get everything past Micro
			# Join everything
			diagnostic = " ".join(data_cells).strip()

	for k in abbreviated_categories:
		diagnostic = re.sub(f"{k}", " M ", diagnostic)

	# Re-stylise anatomies as 'A'
	symbol_list = list(re.findall(r"\w+|[^\w\s]", diagnostic))
	symbol_list = [("A" if i.isnumeric() else i) for i in symbol_list]
	diagnostic = " ".join(symbol_list).strip()

	# Re-stylise ¬ Γ as Γ (nested symbols)
	if re.findall(r"[¬ ]*Γ", diagnostic):
		diagnostic = re.sub(r"[¬ ]*Γ", " Γ", diagnostic)

	drop_indxs = []

	full_row_data = re.findall(r"{(.*)}", value)[0]
	column_data_as_list = re.split(r",", full_row_data)

	symbols = re.findall(r"\w+|[^\w\s]", diagnostic)
	res = [i[0] for i in groupby(symbols)]
	r = " ".join(res)
	vlu = r.strip()

	diagn_as_list = re.split(r"\.", vlu)
	diagn_as_list = list(filter(str.strip, diagn_as_list))

	diagn_set = []

	for item in diagn_as_list:
		if item not in diagn_set:
			if "M" in item:
				diagn_set.append(item)

	diagnostic = ".".join(diagn_set).strip()

	if diagnostic:
		condensed_key_subset.append(key)
		condensed_data[key] = diagnostic

# %%
# Re-number keys in condensed data
df_condensed = df.loc[condensed_key_subset]

data1 = {}
new_keys = list(range(df_condensed.shape[0]))

for new_key, old_key in zip(new_keys, condensed_data.keys()):
	data1[new_key] = condensed_data[old_key]

condensed_data = data1

df_condensed.index = list(condensed_data.keys())

# %%
# Export dictionary as a txt
with open(data_folder / "Processed/condensed_m.txt", "w", encoding="utf-8") as f:
	for key, value in condensed_data.items():
		if value:
			f.write("%s: %s\n" % (int(key), value))
f.close()

# Export dataframe
file_to_save = data_folder / "Processed/condensed_m.csv"
df_condensed.to_csv(file_to_save, index=False)

# %%
