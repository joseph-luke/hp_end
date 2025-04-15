# %%
import time

start = time.time()
# %%
import pandas as pd
import re
from itertools import groupby
from pathlib import Path

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

# Load ROI df
file_to_open_1 = data_folder / "AnatomyPrimes/Anatomy_Prime_Association.xlsx"
df_roi = pd.read_excel(file_to_open_1)

# %%
# Malignant records:
file_to_open_1 = data_folder / "Processed/m_records.csv"
df = pd.read_csv(file_to_open_1)

# %%
# Database (to help subset from data):
file_to_open_1 = data_folder / "Processed/text_corrected_hp.csv"
df2 = pd.read_csv(file_to_open_1, low_memory=False)

desired_records = list(df.SpecimenID)

desired_indices = list(df2.loc[df2.SpecimenID.isin(desired_records)].index)

data2 = {}

for index in desired_indices:
	data2[index] = data[index]

data = data2

# %%
# Re-index malignant data
data1 = {}
new_keys = list(range(df.shape[0]))

for new_key, old_key in zip(new_keys, data.keys()):
	data1[new_key] = data[old_key]

data = data1

df.index = list(data.keys())

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

# %%
# Remove records that do not report non-directional anatomies
no_anatomies = []

for key, value in data.items():
	anatomies = re.findall(r"[0-9]+", value)
	anatomies = list(set(anatomies))
	if len(anatomies) > 0:
		for anatomy in anatomies:
			if (
				int(anatomy) > 48
			):  # 47 was highest directional prime, 577 is directional
				if (
					int(anatomy) != 577
				):  # 47 was highest directional prime, 577 is directional
					pass
	else:
		no_anatomies.append(key)

# %%
for i in no_anatomies:
	del data[i]

# Re-index keys in data
df = df.loc[data.keys()]

data1 = {}
new_keys = list(range(df.shape[0]))

for new_key, old_key in zip(new_keys, data.keys()):
	data1[new_key] = data[old_key]

data = data1

df.index = list(data.keys())

# %%
# Delete encapsulated maps in the row -> retain the more specific map
# Find en(cap)sulated anatomies
A_ = list(df_roi.ROI.str.lower())
A_ = [i.strip() for i in A_]
A_2 = A_.copy()

anaprime = list(df_roi.AnatomicPrime)
anaprime_set = list(sorted(set(df_roi.AnatomicPrime)))

nests = {}

for a1 in A_:
	w1 = anaprime[A_.index(a1)]
	for a2 in A_2:
		w2 = anaprime[A_2.index(a2)]
		if a1 != a2:
			if re.search(a1, a2): # a2 encapsulates a1
				if re.search(a1, a2).span()[0] == 0:
					if anaprime[A_2.index(a2)] != anaprime[A_.index(a1)]:
						if f"{w1} {w2}" not in nests:
							nests[(r"(?<!\d)" + f"{w1} {w2}" + r"(?!\d)")] = f"{str(w2)}"
				elif re.search(a1, a2).span()[0] != 0:
					if anaprime[A_2.index(a2)] != anaprime[A_.index(a1)]:
						if f"{w2} {w1}" not in nests:
							nests[(r"(?<!\d)" + f"{w2} {w1}" + r"(?!\d)")] = f"{str(w2)}"


# %%
denested_data = {}

for key, record in data.items():
	for nest, denest in nests.items():
		if re.findall(nest, record):
			record = re.sub(nest, denest, record)
		denested_data[key] = record

# data2 = data.copy() # just in case we need it
data = denested_data

# %%
joined_data = {}

for key, value in data.items():
	data_columns = []
	if re.search("ClinicalInfo", value):
		data_columns.append("ClinicalInfo")
	if re.search("SpecimenSite", value):
		data_columns.append("SpecimenSite")
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

	indices_to_keep = []

	if "SpecimenSite" in data_columns:
		indices_to_keep.append(data_columns.index("SpecimenSite"))

	if "SpecimenSite" not in data_columns:
		if "ClinicalInfo" in data_columns:
			indices_to_keep.append(data_columns.index("ClinicalInfo"))

	if "Diagnosis" in data_columns:
		indices_to_keep.append(data_columns.index("Diagnosis"))
		if "Comments" in data_columns:
			indices_to_keep.append(data_columns.index("Comments"))

	if "Diagnosis" not in data_columns:
		if "Comments" in data_columns:
			indices_to_keep.append(data_columns.index("Comments"))
		if "Comments" not in data_columns:
			indices_to_keep.append(data_columns.index("MicroDescr"))

	data_columns = [data_columns[i] for i in indices_to_keep]
	data_cells = [data_cells[i] for i in indices_to_keep]

	record = " | ".join(data_cells).strip()

	# Re-stylise ¬ Γ as Γ (nested symbols)
	if re.findall(r"[¬ ]*Γ", record):
		record = re.sub(r"[¬ ]*Γ", " Γ", record)

	# Re-stylise Cs Sa as Cs (nested symbols)
	if re.findall(r"Cs Sa", record):
		record = re.sub(r"Cs Sa", "Cs", record)

	# Delete consecutive occurrences of symbols
	symbols = re.findall(r"\w+|[^\w\s]", record)
	res = [i[0] for i in groupby(symbols)]
	r = " ".join(res)
	vlu = r.strip()

	joined_data[key] = vlu

# %%
all_v_list = []

all_v_joined = []

n_vec_list = []
n_vec = 0

key_vec_list = []

expanded_anatomy_vec_list = []

mapped_anatomies = []
mapped_malignancies = []
mapped_anatomy_primes = []

for key, value in joined_data.items():
	symbol_list = list(re.findall(r"\w+|[^\w\s]", value))
	res = [i[0] for i in groupby(symbol_list)]  # Eliminate consecutive symbols
	symbol_list = res

	anaprimes = list(sorted(set([int(i) for i in symbol_list if i.isnumeric()])))

	tree_node_A1 = []
	tree_node_Anatomy1 = []

	v_anatomies = []
	v_anatomy_primes = []

	for anaprime in anaprimes:
		if anaprime != 577 and anaprime > 48:
			variable_A = str(anaprime)
			partial_A_list1 = [("A" if i == variable_A else i) for i in symbol_list]
			partial_A_list_anatomy1 = [
				(
					list(df_roi.loc[df_roi.AnatomicPrime == int(j)].ROI)[0]
					if j == variable_A
					else j
				)
				for j in symbol_list
			]
			partial_A_list1 = [("" if i.isnumeric() else i) for i in partial_A_list1]
			partial_A_list_anatomy1 = [
				("" if i.isnumeric() and int(i) != 577 and int(i) > 48 else i)
				for i in partial_A_list_anatomy1
			]
			partial_A_list_anatomy1 = list(filter(None, partial_A_list_anatomy1))

			tree_node_A1.append(partial_A_list1)
			tree_node_Anatomy1.append(partial_A_list_anatomy1)

			v_anatomy = list(df_roi.loc[df_roi.AnatomicPrime == int(variable_A)].ROI)[0]
			v_anatomies.append(v_anatomy)

			v_anatomy_primes.append(variable_A)

	for partial_A_list1, partial_A_list_anatomy1 in zip(
		tree_node_A1, tree_node_Anatomy1
	):
		tree_node_A_2 = []
		tree_node_Anatomy_2 = []

		directional_anaprimes = []

		for anaprime in anaprimes:
			if anaprime == 577 or anaprime < 48:
				pass

		if len(directional_anaprimes) > 0:
			partial_A_list2 = [
				("" if i.isnumeric() and int(i) in directional_anaprimes else i)
				for i in partial_A_list1
			]
			partial_A_list2 = list(filter(None, partial_A_list2))

			partial_A_list_anatomy2 = [
				("" if i.isnumeric() else i) for i in partial_A_list_anatomy1
			]
			partial_A_list_anatomy2 = list(filter(None, partial_A_list_anatomy2))

			tree_node_A_2.append(partial_A_list2)
			tree_node_Anatomy_2.append(partial_A_list_anatomy2)

			for anaprime in directional_anaprimes:
				variable_A = str(anaprime)
				partial_A_list2 = [
					("" if i == variable_A else i) for i in partial_A_list1
				]
				partial_A_list_anatomy2 = [
					(
						list(df_roi.loc[df_roi.AnatomicPrime == int(i)].ROI)[0]
						if i == variable_A
						else i
					)
					for i in partial_A_list_anatomy1
				]
				partial_A_list_anatomy2 = [
					("" if i != variable_A and i.isnumeric() else i)
					for i in partial_A_list_anatomy2
				]
				partial_A_list_anatomy2 = list(filter(None, partial_A_list_anatomy2))

				tree_node_A_2.append(partial_A_list2)
				tree_node_Anatomy_2.append(partial_A_list_anatomy2)

		elif len(directional_anaprimes) == 0:
			partial_A_list2 = partial_A_list1
			partial_A_list_anatomy2 = partial_A_list_anatomy1

			tree_node_A_2.append(partial_A_list2)
			tree_node_Anatomy_2.append(partial_A_list_anatomy2)

		# Take derivative wrt M
		v_malignancies = []

		for k in abbreviated_categories:
			m = re.findall(f"{k}", value)
			if m:
				m = m[0]
				if m not in v_malignancies:
					v_malignancies.append(m)

		for malignancy in v_malignancies:
			for partial_A_list2, partial_A_list_anatomy2 in zip(
				tree_node_A_2, tree_node_Anatomy_2
			):
				variable_A = malignancy
				partial_A_list3 = [
					("M" if i == variable_A else i) for i in partial_A_list2
				]
				partial_A_list3 = [
					("m" if i != variable_A and i in v_malignancies else i)
					for i in partial_A_list3
				]

				partial_A_list_anatomy3 = [
					(variable_A if i == variable_A else i)
					for i in partial_A_list_anatomy2
				]
				partial_A_list_anatomy3 = [
					("" if i != variable_A and i in v_malignancies else i)
					for i in partial_A_list_anatomy3
				]
				partial_A_list_anatomy3 = list(filter(None, partial_A_list_anatomy3))

				joined_partial_A_ = " ".join(partial_A_list3).strip()
				joined_partial_A_anatomy = " ".join(partial_A_list_anatomy3).strip()

				symbols_A = re.findall(r"\w+|[^\w\s]", joined_partial_A_)
				symbols_Anatomy = re.findall(r"\w+|[^\w\s]", joined_partial_A_anatomy)

				partial_A_list3 = [i[0] for i in groupby(symbols_A)]
				partial_A_list_anatomy3 = [i[0] for i in groupby(symbols_Anatomy)]

				joined_partial_A_ = " ".join(partial_A_list3).strip()
				joined_partial_A_anatomy = " ".join(partial_A_list_anatomy3).strip()

				all_v_list.append(partial_A_list3)
				all_v_joined.append(joined_partial_A_)
				expanded_anatomy_vec_list.append(joined_partial_A_anatomy)
				n_vec_list.append(n_vec)
				key_vec_list.append(key)

				mapped_malignancies.append(malignancy)

	for a_v in v_anatomies:
		for j in range(len(directional_anaprimes) + 1):
			for k in range(len(v_malignancies)):
				mapped_anatomies.append(a_v)

				n_vec += 1

	for a_v_p in v_anatomy_primes:
		for j in range(len(directional_anaprimes) + 1):
			for k in range(len(v_malignancies)):
				mapped_anatomy_primes.append(a_v_p)

# %%
new_partials = []

for v in all_v_joined:
	split_v = v.split("|")

	score = 0
	m_score = 0

	drop_indices = []

	new_split_v = []

	for col_s in split_v[::-1]:
		if re.findall("A", col_s) and re.findall("M", col_s):
			score += 1
			m_score += 1
			split_y = col_s.split(r"\.")

			board = 0
			drop_indices2 = []

			for sent in split_y[::-1]:

				if re.findall("A", sent) and re.findall("M", sent):
					board += 1

				else:
					board += 0
					if board > 0:
						drop_indices2.append(split_y.index(sent))

			for j in drop_indices2:
				del split_y[j]

			split_y = [m.strip() for m in split_y]

			col_t = r" \. ".join(split_y)

			if score == 1:
				new_split_v.append(col_t)

		else:
			score += 0
			if score > 0:
				pass
			if score == 0:
				if m_score == 0:
					if re.findall("A", col_s) or re.findall("M", col_s):
						m_score += 1
						new_split_v.append(col_s)
				if m_score > 0:
					if re.findall("A", col_s) or re.findall("M", col_s):
						new_split_v.append(col_s)

	new_split_v.reverse()

	split_v = new_split_v

	split_v = [i.strip() for i in split_v]

	joined_partial_A_ = " | ".join(split_v).strip()

	new_partials.append(joined_partial_A_)


all_v_joined = new_partials


# Remove repeated columns
new_partials = []

for v in all_v_joined:
	split_v = v.split("|")
	split_v = [i.strip() for i in split_v]
	split_v = [i[0] for i in groupby(split_v)]

	split_v = [i.strip() for i in split_v]

	joined_partial_A_ = " | ".join(split_v).strip()

	new_partials.append(joined_partial_A_)


all_v_joined = new_partials

# %%
# Export dictionary as a txt
with open(data_folder / "Processed/condensed_am.txt", "w", encoding="utf-8") as f:
	for key, value in enumerate(all_v_joined):
		f.write("%s: %s\n" % (int(key), value))
f.close()

# Export dataframe
file_to_save = data_folder / "Processed/condensed_am.csv"
df.to_csv(file_to_save, index=False)

# %%
# Export for future A-M record localization
potential_maps_df = pd.DataFrame(
	{
		"record": key_vec_list,
		"anatomy": mapped_anatomies,
		"anatomyprime": mapped_anatomy_primes,
		"malignancy": mapped_malignancies,
	}
)

file_to_save = data_folder / "Processed/potential_am.csv"
potential_maps_df.to_csv(file_to_save, index=False)


# %%
end = time.time()
print("The program took: ", (end - start) / 60, "min", "to run.")

# %%
