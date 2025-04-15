# %%
import time

start = time.time()

# %%
import numpy as np
import ahocorasick
import pandas as pd
from bisect import bisect_left
from pathlib import Path


def take_closest(sorted_list, number):
	pos = bisect_left(sorted_list, number)
	if pos == 0:
		return sorted_list[0]
	if pos == len(sorted_list):
		return sorted_list[-1]
	after = sorted_list[pos]
	return after


def distance_to_closest(sorted_list, number):
	pos = bisect_left(sorted_list, number)
	if pos == 0:
		return sorted_list[0] - number
	if pos == len(sorted_list):
		return sorted_list[-1] - number
	after = sorted_list[pos]
	return after - number


automaton = ahocorasick.Automaton()
automaton2 = ahocorasick.Automaton()
automaton3 = ahocorasick.Automaton()
automaton4 = ahocorasick.Automaton()
automaton5 = ahocorasick.Automaton()
automaton6 = ahocorasick.Automaton()
automaton7 = ahocorasick.Automaton()
automaton8 = ahocorasick.Automaton()
automaton9 = ahocorasick.Automaton()
automaton10 = ahocorasick.Automaton()
automaton11 = ahocorasick.Automaton()
automaton12 = ahocorasick.Automaton()
automaton13 = ahocorasick.Automaton()
automaton14 = ahocorasick.Automaton()
automaton15 = ahocorasick.Automaton()
automaton16 = ahocorasick.Automaton()
automaton17 = ahocorasick.Automaton()

# %%
data_folder = Path(
	"C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"
)

file_to_open_1 = data_folder / "Processed/text_corrected_hp.csv"
df = pd.read_csv(file_to_open_1, low_memory=False)
df = df.iloc[0:10]

df_ticked = df.copy(deep=True)

for i in range(1, (df_ticked.shape[1] * 2) + 1, 2):
	df_ticked.insert(i, "coord_tick", ["coord_tick"] * df_ticked.shape[0], allow_duplicates=True)

# %%
data = df_ticked.to_csv(
	None, ",", na_rep=None, index=False, index_label=None, lineterminator=None
)

# %%
categories = (
	"carcinosarcoma",
	"carcinoma",
	"sarcoma",
	"lymphoma",
	"leukemia",
	"blastoma",
	"glioma",
	"melanoma",
	"mesothelioma",
	"teratoma",
	"seminoma",
	"astrocytoma",
	"peripheral nerve sheath tumor",
	"myeloma",
)

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
]

# %%
full_data_indices = []
full_data_indices = np.array(full_data_indices)

wanted_cell_data_indices_lb = []
wanted_cell_data_indices_lb = np.array(wanted_cell_data_indices_lb)


wanted_cell_data_indices_ub = []
wanted_cell_data_indices_ub = np.array(wanted_cell_data_indices_ub)

cat_data_indices_dict = {}

full_dict = {}
full_row_indices = []

full_col_indices = []
full_col_indices = np.array(full_col_indices)

for (category), (ctg) in zip(categories, abbreviated_categories):
	cat = category

	automaton = ahocorasick.Automaton()
	automaton.add_word(cat, (1, cat))

	automaton.make_automaton()

	automaton2.add_word(",coord_tick", (1, ",coord_tick"))
	automaton2.make_automaton()

	# These are all coordinate data indices
	coord_ticks = np.array([i[0] for i in list(automaton2.iter_long(data))])

	# These are the data indices where the terms we want start
	data_indices = np.array([i[0] for i in list(automaton.iter_long(data))])

	# These are the coord_ticks closest to the terms
	wanted_coord_ticks = np.array(
		list(sorted(set(take_closest(coord_ticks, i) for i in data_indices)))
	)

	# These are the indices of the coord_ticks, which will tell us the column
	# the term came from
	col_indices = np.where(np.in1d(coord_ticks, wanted_coord_ticks))

	full_col_indices = np.concatenate((full_col_indices, col_indices[0])).astype(int)

	col_indices = col_indices[0]
	col_indicesb = col_indices - 1

	# This creates a dictionary for each column that contains the indices for the col
	column_map = {}

	for i in range(0, df.shape[1]):
		column_indices = [i]
		for m in range(1, df.shape[0] + 1):
			column_indices.append(m * df.shape[1] + i)
		column_map[df.columns[i]] = np.array(column_indices)

	# This tells us which column the malignant term reference came from
	column_with_diagnosis = []
	for w in col_indices:
		for key in column_map:
			if w in column_map[key]:
				column_with_diagnosis.append(key)
				break  # Break immediately we find the column

	# This tells us which row the malignant term reference came from
	row_ends = column_map[df.columns[df.shape[1] - 1]]
	row_ticks = np.array(list(sorted(take_closest(row_ends, i) for i in col_indices)))

	# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
	row_indices = []
	for i in row_ticks:
		r = np.where(row_ends == i)[0][0]
		row_indices.append(r)

	row_indices = np.array(row_indices) - 1

	list_row_indices = list(row_indices)

	full_row_indices += list_row_indices

	len(row_indices)  # Missing 8 spellings -> correct spelling/capitalisation!

	# Create a dictionary that summarises the diagnoses rows and columns
	diagnosis_location_dict = {}
	for i in sorted(set(row_indices)):
		diagnosis_location_dict[i] = []

	for i in range(len(row_indices)):
		diagnosis_and_loc = column_with_diagnosis[i] + "-" + str(data_indices[i])
		diagnosis_location_dict[row_indices[i]].append(diagnosis_and_loc)

	full_dict[category] = diagnosis_location_dict

	# These are the indices of the coord_ticks, which will tell us the column
	# the term came from
	col_indices = np.where(np.in1d(coord_ticks, wanted_coord_ticks))
	col_indices = col_indices[0]  # far bound
	col_indicesb = col_indices - 1  # near bound

	# These give me all the malignant cell boundaries (and I will parse through)
	wanted_cell_data_indices_near_b = coord_ticks[col_indicesb]
	wanted_cell_data_indices_far_b = coord_ticks[col_indices]

	wanted_cell_data_indices_lb = np.concatenate(
		(wanted_cell_data_indices_lb, wanted_cell_data_indices_near_b)
	)
	wanted_cell_data_indices_ub = np.concatenate(
		(wanted_cell_data_indices_ub, wanted_cell_data_indices_far_b)
	)

	full_data_indices = np.concatenate((full_data_indices, data_indices))

	cat_data_indices_dict[ctg] = data_indices

full_cols = []
for w in full_col_indices:
	for key in column_map:
		if w in column_map[key]:
			full_cols.append(key)

# %%
# Stop terms
alphabet = [
	"a",
	"b",
	"c",
	"d",
	"e",
	"f",
	"g",
	"h",
	"i",
	"j",
	"k",
	"l",
	"m",
	"n",
	"o",
	"p",
	"q",
	"r",
	"s",
	"t",
	"u",
	"v",
	"w",
	"x",
	"y",
	"z",
]

alpha_bullet = []
alpha_bullet += ["\n" + i + "]" for i in alphabet]
alpha_bullet += ["\n" + i + ")" for i in alphabet]

alpha_num = []
alpha_num += ["\n" + str(i) + "]" for i in range(0, 101)]
alpha_num += ["\n" + str(i) + ")" for i in range(0, 101)]

alpha_bullet += alpha_num

alpha_bullet += [".", ";", "\n-"]


for index, word in enumerate(alpha_bullet):
	automaton3.add_word(word, (index, word))

automaton3.make_automaton()

# These are the indices for the stop markers
stop_indices = np.array([i[0] for i in list(automaton3.iter_long(data))])

# %%
# These are the coord_ticks closest to the terms
wanted_stop_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in stop_indices)))
)

col_stop_indices = np.where(np.in1d(coord_ticks, wanted_stop_coord_ticks))

# This tells us which column the stop marker came from
column_with_stop_marker = []
for w in col_stop_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_stop_marker.append(key)
			break  # Break immediately we find the column

# This tells us which row the stop marker came from
row_stop_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_stop_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_stop_indices = []
for i in row_stop_ticks:
	r = np.where(row_ends == i)[0][0]
	row_stop_indices.append(r)

row_stop_indices = np.array(row_stop_indices) - 1

# Create a dictionary that summarises the stop marker rows and columns
stop_marker_location_dict = {}
for i in sorted(set(row_stop_indices)):
	stop_marker_location_dict[i] = []

for i in range(len(row_stop_indices)):
	stop_and_loc = column_with_stop_marker[i] + "-" + str(stop_indices[i])
	stop_marker_location_dict[row_stop_indices[i]].append(stop_and_loc)

# %%
# Global negators
global_negators = [
	" no malignancy",
	" no evidence of malignancy",
	" not malignant",
	" negative for malignancy",
	"not seen",
	"not found",
	"not identified",
	"not observed",
	"is absent",
	"was absent",
	"are absent",
	"be absent",
	"no evidence of",
	"rebiopsy",
	"re-biopsy",
	"re - biopsy" "\nnone",
	" none",
]

for index, word in enumerate(global_negators):
	automaton4.add_word(word, (index, word))
automaton4.make_automaton()

# These are the indices for the global_negators
global_negator_indices = np.array(
	[i[0] for i in list(automaton4.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_global_negator_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in global_negator_indices)))
)

col_global_negator_indices = np.where(
	np.in1d(coord_ticks, wanted_global_negator_coord_ticks)
)

# This tells us which column the global_negator came from
column_with_global_negator = []
for w in col_global_negator_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_global_negator.append(key)
			break  # Break immediately we find the column

# This tells us which row the global_negator came from
row_global_negator_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_global_negator_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_global_negator_indices = []
for i in row_global_negator_ticks:
	r = np.where(row_ends == i)[0][0]
	row_global_negator_indices.append(r)

row_global_negator_indices = np.array(row_global_negator_indices) - 1

# Create a dictionary that summarises the global_negator rows and columns
global_negator_location_dict = {}
for i in sorted(set(row_global_negator_indices)):
	global_negator_location_dict[i] = []

for i in range(len(row_global_negator_indices)):
	global_negator_location_dict[row_global_negator_indices[i]].append(
		column_with_global_negator[i]
	)

# This will give us the cell limits that contain global negators
col_global_negator_indices = col_global_negator_indices[0]  # far bound
col_global_negator_indicesb = col_global_negator_indices - 1  # near bound

# These give me all the malignant cell boundaries (and I will parse through)
wanted_cell_glob_neg_indices_near_b = coord_ticks[col_global_negator_indicesb]
wanted_cell_glob_neg_indices_far_b = coord_ticks[col_global_negator_indices]


# %%
# Specific negators
precise_negators = [
	"negative",
	" not ",
	" no ",
	"-no ",
	";no ",
	"\nnot",
	"-not ",
	"nno ",
	"\nno ",
	"nno ",
	" without ",
	"no ",
	"\nwithout",
	"free of",
	"precluding",
	"preclude",
	"excluding",
	"exclude",
	"negate",
	"\nnegate",
	"insufficient",
]

for index, word in enumerate(precise_negators):
	automaton5.add_word(word, (index, word))
automaton5.make_automaton()

# These are the indices for the precise_negators
precise_negator_indices = np.array(
	[i[0] for i in list(automaton5.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_precise_negator_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in precise_negator_indices)))
)

col_precise_negator_indices = np.where(
	np.in1d(coord_ticks, wanted_precise_negator_coord_ticks)
)

# This tells us which column the precise_negator came from
column_with_precise_negator = []
for w in col_precise_negator_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_precise_negator.append(key)
			break  # Break immediately we find the column

# This tells us which row the precise_negator came from
row_precise_negator_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_precise_negator_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_precise_negator_indices = []
for i in row_precise_negator_ticks:
	r = np.where(row_ends == i)[0][0]
	row_precise_negator_indices.append(r)

row_precise_negator_indices = np.array(row_precise_negator_indices) - 1

# Create a dictionary that summarises the precise_negator rows and columns
precise_negator_location_dict = {}
for i in sorted(set(row_precise_negator_indices)):
	precise_negator_location_dict[i] = []

for i in range(len(row_precise_negator_indices)):
	precise_negator_location_dict[row_precise_negator_indices[i]].append(
		column_with_precise_negator[i] + ": ¬"
	)

# %%
# Specific posators
precise_posators = ["however", "but ", "nevertheless", "nbut", "although"]

for index, word in enumerate(precise_posators):
	automaton8.add_word(word, (index, word))
automaton8.make_automaton()

# These are the indices for the precise_posators
precise_posator_indices = np.array(
	[i[0] for i in list(automaton8.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_precise_posator_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in precise_posator_indices)))
)

col_precise_posator_indices = np.where(
	np.in1d(coord_ticks, wanted_precise_posator_coord_ticks)
)

# This tells us which column the precise_posator came from
column_with_precise_posator = []
for w in col_precise_posator_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_precise_posator.append(key)
			break  # Break immediately we find the column

# This tells us which row the precise_posator came from
row_precise_posator_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_precise_posator_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_precise_posator_indices = []
for i in row_precise_posator_ticks:
	r = np.where(row_ends == i)[0][0]
	row_precise_posator_indices.append(r)

row_precise_posator_indices = np.array(row_precise_posator_indices) - 1

# Create a dictionary that summarises the precise_posator rows and columns
precise_posator_location_dict = {}
for i in sorted(set(row_precise_posator_indices)):
	precise_posator_location_dict[i] = []

for i in range(len(row_precise_posator_indices)):
	precise_posator_location_dict[row_precise_posator_indices[i]].append(
		column_with_precise_posator[i] + ": ¬"
	)


# %%
# Inverse terms
inverse_terms = ["benign", "normal"]

for index, word in enumerate(inverse_terms):
	automaton6.add_word(word, (index, word))
automaton6.make_automaton()

# These are the indices for the inverse_terms
inverse_term_indices = np.array([i[0] for i in list(automaton6.iter_long(data))])

# These are the coord_ticks closest to the terms
wanted_inverse_term_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in inverse_term_indices)))
)

col_inverse_term_indices = np.where(
	np.in1d(coord_ticks, wanted_inverse_term_coord_ticks)
)

# This tells us which column the inverse_term came from
column_with_inverse_term = []
for w in col_inverse_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_inverse_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the inverse_term came from
row_inverse_term_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_inverse_term_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_inverse_term_indices = []
for i in row_inverse_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_inverse_term_indices.append(r)

row_inverse_term_indices = np.array(row_inverse_term_indices) - 1

# Create a dictionary that summarises the inverse_term rows and columns
inverse_term_location_dict = {}
for i in sorted(set(row_inverse_term_indices)):
	inverse_term_location_dict[i] = []

for i in range(len(row_inverse_term_indices)):
	inverse_term_location_dict[row_inverse_term_indices[i]].append(
		column_with_inverse_term[i]
	)

# %%
# Margin terms
margin_term = "margin"

automaton9.add_word(margin_term, (1, margin_term))

automaton9.make_automaton()

# These are the indices for the margin_terms
margin_term_indices = np.array([i[0] for i in list(automaton9.iter_long(data))])

# These are the coord_ticks closest to the terms
wanted_margin_term_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in margin_term_indices)))
)

col_margin_term_indices = np.where(np.in1d(coord_ticks, wanted_margin_term_coord_ticks))

# This tells us which column the margin_term came from
column_with_margin_term = []
for w in col_margin_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_margin_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the margin_term came from
row_margin_term_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_margin_term_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_margin_term_indices = []
for i in row_margin_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_margin_term_indices.append(r)

row_margin_term_indices = np.array(row_margin_term_indices) - 1

# Create a dictionary that summarises the margin_term rows and columns
margin_term_location_dict = {}
for i in sorted(set(row_margin_term_indices)):
	margin_term_location_dict[i] = []

for i in range(len(row_margin_term_indices)):
	margin_term_location_dict[row_margin_term_indices[i]].append(
		column_with_margin_term[i]
	)

# %%
# Cancer-anatomy mapping terms
mapping_terms = [" is in ", " show ", " of the ", ":", "- "]

for index, word in enumerate(mapping_terms):
	automaton11.add_word(word, (index, word))
automaton11.make_automaton()

# These are the indices for the mapping_terms
mapping_term_indices = np.array(
	[i[0] for i in list(automaton11.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_mapping_term_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in mapping_term_indices)))
)

col_mapping_term_indices = np.where(
	np.in1d(coord_ticks, wanted_mapping_term_coord_ticks)
)

# This tells us which column the mapping_term came from
column_with_mapping_term = []
for w in col_mapping_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_mapping_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the mapping_term came from
row_mapping_term_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_mapping_term_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_mapping_term_indices = []
for i in row_mapping_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_mapping_term_indices.append(r)

row_mapping_term_indices = np.array(row_mapping_term_indices) - 1

# Create a dictionary that summarises the mapping_term rows and columns
mapping_term_location_dict = {}
for i in sorted(set(row_mapping_term_indices)):
	mapping_term_location_dict[i] = []

for i in range(len(row_mapping_term_indices)):
	mapping_term_location_dict[row_mapping_term_indices[i]].append(
		column_with_mapping_term[i]
	)

# %%
# Ambivalence terms
ambivalence_terms = [" or ", "\nor ", " or ", " vs ", "\nvs ", "versus"]

for index, word in enumerate(ambivalence_terms):
	automaton13.add_word(word, (index, word))
automaton13.make_automaton()

# These are the indices for the ambivalence_terms
ambivalence_term_indices = np.array(
	[i[0] for i in list(automaton13.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_ambivalence_term_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in ambivalence_term_indices)))
)

col_ambivalence_term_indices = np.where(
	np.in1d(coord_ticks, wanted_ambivalence_term_coord_ticks)
)

# This tells us which column the ambivalence_term came from
column_with_ambivalence_term = []
for w in col_ambivalence_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_ambivalence_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the ambivalence_term came from
row_ambivalence_term_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_ambivalence_term_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_ambivalence_term_indices = []
for i in row_ambivalence_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_ambivalence_term_indices.append(r)

row_ambivalence_term_indices = np.array(row_ambivalence_term_indices) - 1

# Create a dictionary that summarises the ambivalence_term rows and columns
ambivalence_term_location_dict = {}
for i in sorted(set(row_ambivalence_term_indices)):
	ambivalence_term_location_dict[i] = []

for i in range(len(row_ambivalence_term_indices)):
	ambivalence_term_location_dict[row_ambivalence_term_indices[i]].append(
		column_with_ambivalence_term[i]
	)

# %%
# Specific benign-indicator terms
benign_indicator_terms = ["fibroadenocarcinoma", "mature"]

for index, word in enumerate(benign_indicator_terms):
	automaton12.add_word(word, (index, word))
automaton12.make_automaton()

# These are the indices for the benign_indicator_terms
benign_indicator_term_indices = np.array(
	[i[0] for i in list(automaton12.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_benign_indicator_term_coord_ticks = np.array(
	list(
		sorted(set(take_closest(coord_ticks, i) for i in benign_indicator_term_indices))
	)
)

col_benign_indicator_term_indices = np.where(
	np.in1d(coord_ticks, wanted_benign_indicator_term_coord_ticks)
)

# This tells us which column the benign_indicator_term came from
column_with_benign_indicator_term = []
for w in col_benign_indicator_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_benign_indicator_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the benign_indicator_term came from
row_benign_indicator_term_ticks = np.array(
	list(
		sorted(take_closest(row_ends, i) for i in col_benign_indicator_term_indices[0])
	)
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_benign_indicator_term_indices = []
for i in row_benign_indicator_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_benign_indicator_term_indices.append(r)

row_benign_indicator_term_indices = np.array(row_benign_indicator_term_indices) - 1

# Create a dictionary that summarises the benign_indicator_term rows and columns
benign_indicator_term_location_dict = {}
for i in sorted(set(row_benign_indicator_term_indices)):
	benign_indicator_term_location_dict[i] = []

for i in range(len(row_benign_indicator_term_indices)):
	benign_indicator_term_location_dict[row_benign_indicator_term_indices[i]].append(
		column_with_benign_indicator_term[i]
	)


# %%
# Specific premalignant-indicator terms
premalignant_indicator_terms = [
	"in situ",
	"metaplasia",
	"pre-malignant",
	"premalignant",
	"pre malignant",
	"pre-malignancy",
	"premalignancy",
	"pre malignancy",
	"intraepithelial",
	"non-infiltrating",
	"noninfiltrating",
	"non-invasive",
	"noninvasive",
	"pre-cancerous",
	"precancerous",
	"pre cancerous",
	"in-situ",
	"insitu",
	"in - situ",
	"in -situ",
]

for index, word in enumerate(premalignant_indicator_terms):
	automaton17.add_word(word, (index, word))
automaton17.make_automaton()

# These are the indices for the premalignant_indicator_terms
premalignant_indicator_term_indices = np.array(
	[i[0] for i in list(automaton17.iter(data, ignore_white_space=False))]
)

# These are the coord_ticks closest to the terms
wanted_premalignant_indicator_term_coord_ticks = np.array(
	list(
		sorted(
			set(
				take_closest(coord_ticks, i)
				for i in premalignant_indicator_term_indices
			)
		)
	)
)

col_premalignant_indicator_term_indices = np.where(
	np.in1d(coord_ticks, wanted_premalignant_indicator_term_coord_ticks)
)

# This tells us which column the premalignant_indicator_term came from
column_with_premalignant_indicator_term = []
for w in col_premalignant_indicator_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_premalignant_indicator_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the premalignant_indicator_term came from
row_premalignant_indicator_term_ticks = np.array(
	list(
		sorted(
			take_closest(row_ends, i)
			for i in col_premalignant_indicator_term_indices[0]
		)
	)
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_premalignant_indicator_term_indices = []
for i in row_premalignant_indicator_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_premalignant_indicator_term_indices.append(r)

row_premalignant_indicator_term_indices = (
	np.array(row_premalignant_indicator_term_indices) - 1
)

# Create a dictionary that summarises the premalignant_indicator_term rows and columns
premalignant_indicator_term_location_dict = {}
for i in sorted(set(row_premalignant_indicator_term_indices)):
	premalignant_indicator_term_location_dict[i] = []

for i in range(len(row_premalignant_indicator_term_indices)):
	premalignant_indicator_term_location_dict[
		row_premalignant_indicator_term_indices[i]
	].append(column_with_premalignant_indicator_term[i])

# %%
# Differential terms
differential_terms = [
	"differential",
	"ddx",
	"between",
	"possib",
	"maybe",
	"rule out",
	"ruled out",
	"ruling out",
	"r/o",
	" ro ",
	"could",
	"might",
	"consider",
	"perhap",
	"probabl",
	"uncertain",
	"may be",
	"susp"
]

for index, word in enumerate(differential_terms):
	automaton14.add_word(word, (index, word))
automaton14.make_automaton()

# These are the indices for the differential_terms
differential_term_indices = np.array([i[0] for i in list(automaton14.iter_long(data))])

# These are the coord_ticks closest to the terms
wanted_differential_term_coord_ticks = np.array(
	list(sorted(set(take_closest(coord_ticks, i) for i in differential_term_indices)))
)

col_differential_term_indices = np.where(
	np.in1d(coord_ticks, wanted_differential_term_coord_ticks)
)

# This tells us which column the differential_term came from
column_with_differential_term = []
for w in col_differential_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_differential_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the differential_term came from
row_differential_term_ticks = np.array(
	list(sorted(take_closest(row_ends, i) for i in col_differential_term_indices[0]))
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_differential_term_indices = []
for i in row_differential_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_differential_term_indices.append(r)

row_differential_term_indices = np.array(row_differential_term_indices) - 1

# Create a dictionary that summarises the differential_term rows and columns
differential_term_location_dict = {}
for i in sorted(set(row_differential_term_indices)):
	differential_term_location_dict[i] = []

for i in range(len(row_differential_term_indices)):
	differential_term_location_dict[row_differential_term_indices[i]].append(
		column_with_differential_term[i]
	)

# %%
# Malignancy terms
unspecified_malignancy_terms = [
	"malignant neoplas",
	"malignant tumor",
	"round cell malignancy",
	"malignant round cell tumor",
	"malignant small round cell",
	"small round cell malignant",
	"spindle cell malignancy",
	"malignant spindle cell",
	"malignant growth",
	"small blue cell",
	"round blue cell",
	"unspecified malignancy",
	"malignancy unspecified",
	"carcinoid",
	"small cell neoplas",
	"small cell tumor",
	"small cell malignancy",
	"high grade neoplasm",
	"poorly differentiated neoplasm",
	"high grade malignancy",
	"high grade tumor",
	"poorly differentiated malignancy",
	"poorly differentiated tumor",
]

for index, word in enumerate(unspecified_malignancy_terms):
	automaton15.add_word(word, (index, word))
automaton15.make_automaton()

# These are the indices for the unspecified_malignancy_terms
unspecified_malignancy_term_indices = np.array(
	[i[0] for i in list(automaton15.iter_long(data))]
)

# These are the coord_ticks closest to the terms
wanted_unspecified_malignancy_term_coord_ticks = np.array(
	list(
		sorted(
			set(
				take_closest(coord_ticks, i)
				for i in unspecified_malignancy_term_indices
			)
		)
	)
)

col_unspecified_malignancy_term_indices = np.where(
	np.in1d(coord_ticks, wanted_unspecified_malignancy_term_coord_ticks)
)

# This tells us which column the unspecified_malignancy_term came from
column_with_unspecified_malignancy_term = []
for w in col_unspecified_malignancy_term_indices[0]:
	for key in column_map:
		if w in column_map[key]:
			column_with_unspecified_malignancy_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the unspecified_malignancy_term came from
row_unspecified_malignancy_term_ticks = np.array(
	list(
		sorted(
			take_closest(row_ends, i)
			for i in col_unspecified_malignancy_term_indices[0]
		)
	)
)

# Couldn't use np.where(np.in1d) because it removed duplicates for some reason
row_unspecified_malignancy_term_indices = []
for i in row_unspecified_malignancy_term_ticks:
	r = np.where(row_ends == i)[0][0]
	row_unspecified_malignancy_term_indices.append(r)

row_unspecified_malignancy_term_indices = (
	np.array(row_unspecified_malignancy_term_indices) - 1
)

# Create a dictionary that summarises the unspecified_malignancy_term rows and columns
unspecified_malignancy_term_location_dict = {}
for i in sorted(set(row_unspecified_malignancy_term_indices)):
	unspecified_malignancy_term_location_dict[i] = []

for i in range(len(row_unspecified_malignancy_term_indices)):
	unspecified_malignancy_term_location_dict[
		row_unspecified_malignancy_term_indices[i]
	].append(column_with_unspecified_malignancy_term[i])

# This will give us the cell limits that contain unspecified malignancy terms
col_unspecified_malignancy_term_indices = col_unspecified_malignancy_term_indices[
	0
]  # far bound
col_unspecified_malignancy_term_indicesb = (
	col_unspecified_malignancy_term_indices - 1
)  # near bound

# These give me all the unspecified malignancy cell boundaries (and I will parse through)
wanted_cell_unspecified_malignancy_indices_near_b = coord_ticks[
	col_unspecified_malignancy_term_indicesb
]
wanted_cell_unspecified_malignancy_indices_far_b = coord_ticks[
	col_unspecified_malignancy_term_indices
]

# %%
# Load anatomy-prime association df
file_to_open_1 = data_folder / "AnatomyPrimes/Anatomy_Prime_Association.xlsx"
df_roi = pd.read_excel(file_to_open_1)

for column in df_roi.columns:
	try:
		df_roi[column] = df_roi[column].str.lower()
	except:
		pass

for column in df_roi.columns:
	try:
		df_roi[column] = df_roi[column].str.strip()
	except:
		pass

# %%
ROIs = list(df_roi["ROI"])
primes = list(df_roi["AnatomicPrime"])

# %%
full_anapri_data_indices = []
full_anapri_data_indices = np.array(full_anapri_data_indices)

wanted_cell_anapri_data_indices_lb = []
wanted_cell_anapri_data_indices_lb = np.array(wanted_cell_anapri_data_indices_lb)


wanted_cell_anapri_data_indices_ub = []
wanted_cell_anapri_data_indices_ub = np.array(wanted_cell_anapri_data_indices_ub)

roi_anapri_data_indices_dict = {}
for prime in primes:
	roi_anapri_data_indices_dict[prime] = np.array([])

full_anapri_dict = {}
full_anapri_row_indices = []

full_anapri_col_indices = []
full_anapri_col_indices = np.array(full_anapri_col_indices)


for roi, prime in zip(ROIs, primes):
	print(roi)
	roi = roi.strip()

	if len(roi) < 4:
		roi = " " + roi + " "

	if roi == "cardia":
		roi = " cardia"

	if roi == "face":
		roi = " face "

	if roi == "axilla":
		roi = " axilla"

	if roi == "wing":
		roi = " wing"

	if roi == "nose":
		roi = " nose"

	if roi == "oral":
		roi = " oral"

	if roi == "pons":
		roi = " pons"

	if roi == "neural":
		roi = " neural "

	if roi == "lymph":
		roi = " lymph "

	if roi == "colon":
		roi = " colon "

	if roi == "cranial":
		roi = " cranial"

	automaton10 = ahocorasick.Automaton()
	automaton10.add_word(roi, (1, roi))
	automaton10.make_automaton()

	automaton16.add_word(",coord_tick", (1, ",coord_tick"))
	automaton16.make_automaton()

	# These are all coordinate data indices
	coord_ticks = np.array([i[0] for i in list(automaton16.iter_long(data))])

	# These are the data indices where the terms we want start
	anapri_data_indices = np.array([i[0] for i in list(automaton10.iter_long(data))])

	# These are the coord_ticks closest to the terms
	wanted_coord_ticks = np.array(
		list(sorted(set(take_closest(coord_ticks, i) for i in anapri_data_indices)))
	)

	# These are the indices of the coord_ticks, which will tell us the column
	# the term came from
	col_indices = np.where(np.in1d(coord_ticks, wanted_coord_ticks))

	full_anapri_col_indices = np.concatenate(
		(full_anapri_col_indices, col_indices[0])
	).astype(int)

	col_indices = col_indices[0]
	col_indicesb = col_indices - 1

	# This creates a dictionary for each column that contains the indices for the col
	column_map = {}

	for i in range(0, df.shape[1]):
		column_indices = [i]
		for m in range(1, df.shape[0] + 1):
			column_indices.append(m * df.shape[1] + i)
		column_map[df.columns[i]] = np.array(column_indices)

	# This tells us which column the malignant term reference came from
	column_with_diagnosis = []
	for w in col_indices:
		for key in column_map:
			if w in column_map[key]:
				column_with_diagnosis.append(key)
				break  # Break immediately we find the column

	# This tells us which row the malignant term reference came from
	row_ends = column_map[df.columns[df.shape[1] - 1]]
	row_ticks = np.array(list(sorted(take_closest(row_ends, i) for i in col_indices)))

	# Couldn't use np.where(np.in1d) because it removed dupliroies for some reason
	row_indices = []
	for i in row_ticks:
		r = np.where(row_ends == i)[0][0]
		row_indices.append(r)

	row_indices = np.array(row_indices) - 1

	list_row_indices = list(row_indices)

	full_anapri_row_indices += list_row_indices

	len(row_indices)

	# Create a dictionary that summarises the diagnoses rows and columns
	diagnosis_location_dict = {}
	for i in sorted(set(row_indices)):
		diagnosis_location_dict[i] = []

	for i in range(len(row_indices)):
		diagnosis_and_loc = column_with_diagnosis[i] + "-" + str(anapri_data_indices[i])
		diagnosis_location_dict[row_indices[i]].append(diagnosis_and_loc)

	full_anapri_dict[roi] = diagnosis_location_dict

	# These are the indices of the coord_ticks, which will tell us the column
	# the term came from
	col_indices = np.where(np.in1d(coord_ticks, wanted_coord_ticks))
	col_indices = col_indices[0]  # far bound
	col_indicesb = col_indices - 1  # near bound

	# These give me all the malignant cell boundaries (and I will parse through)
	wanted_cell_anapri_data_indices_near_b = coord_ticks[col_indicesb]
	wanted_cell_anapri_data_indices_far_b = coord_ticks[col_indices]

	wanted_cell_anapri_data_indices_lb = np.concatenate(
		(wanted_cell_anapri_data_indices_lb, wanted_cell_anapri_data_indices_near_b)
	)
	wanted_cell_anapri_data_indices_ub = np.concatenate(
		(wanted_cell_anapri_data_indices_ub, wanted_cell_anapri_data_indices_far_b)
	)

	full_anapri_data_indices = np.concatenate(
		(full_anapri_data_indices, anapri_data_indices)
	)

	roi_anapri_data_indices_dict[prime] = np.append(
		roi_anapri_data_indices_dict[prime], anapri_data_indices
	)

	print(len(roi_anapri_data_indices_dict[prime]))

full_anapri_cols = []
for w in full_anapri_col_indices:
	for key in column_map:
		if w in column_map[key]:
			full_anapri_cols.append(key)

# %%
full_data_indices = np.sort(full_data_indices)

# Add cells with global negators to the 'wanted' cell list
missing_indices1 = np.where(
	np.in1d(
		wanted_cell_glob_neg_indices_near_b, wanted_cell_data_indices_lb, invert=True
	)
)

missing_glob_cell_values_lb = wanted_cell_glob_neg_indices_near_b[
	missing_indices1
].astype(int)
missing_glob_cell_values_ub = wanted_cell_glob_neg_indices_far_b[
	missing_indices1
].astype(int)

step1_wanted_cell_data_indices_lb = np.concatenate(
	(wanted_cell_data_indices_lb, missing_glob_cell_values_lb)
)
step1_wanted_cell_data_indices_ub = np.concatenate(
	(wanted_cell_data_indices_ub, missing_glob_cell_values_ub)
)

step1_wanted_cell_data_indices_lb = np.sort(step1_wanted_cell_data_indices_lb)
step1_wanted_cell_data_indices_ub = np.sort(step1_wanted_cell_data_indices_ub)

step1_wanted_cell_data_indices_lb = sorted(set(step1_wanted_cell_data_indices_lb))
step1_wanted_cell_data_indices_ub = sorted(set(step1_wanted_cell_data_indices_ub))

# Add cells with unspecified malignancy terms to the 'wanted' cell list
missing_indices2 = np.where(
	np.in1d(
		wanted_cell_unspecified_malignancy_indices_near_b,
		step1_wanted_cell_data_indices_lb,
		invert=True,
	)
)

missing_gen_m_cell_values_lb = wanted_cell_unspecified_malignancy_indices_near_b[
	missing_indices2
].astype(int)
missing_gen_m_cell_values_ub = wanted_cell_unspecified_malignancy_indices_far_b[
	missing_indices2
].astype(int)

step2_wanted_cell_data_indices_lb = np.concatenate(
	(step1_wanted_cell_data_indices_lb, missing_gen_m_cell_values_lb)
)
step2_wanted_cell_data_indices_ub = np.concatenate(
	(step1_wanted_cell_data_indices_ub, missing_gen_m_cell_values_ub)
)

step2_wanted_cell_data_indices_lb = np.sort(step2_wanted_cell_data_indices_lb)
step2_wanted_cell_data_indices_ub = np.sort(step2_wanted_cell_data_indices_ub)

step2_wanted_cell_data_indices_lb = sorted(set(step2_wanted_cell_data_indices_lb))
step2_wanted_cell_data_indices_ub = sorted(set(step2_wanted_cell_data_indices_ub))

# Add cells with anatomy terms to the 'wanted' cell list
missing_indices3 = np.where(
	np.in1d(
		wanted_cell_anapri_data_indices_lb,
		step2_wanted_cell_data_indices_lb,
		invert=True,
	)
)

missing_anapri_cell_values_lb = wanted_cell_anapri_data_indices_lb[
	missing_indices3
].astype(int)
missing_anapri_cell_values_ub = wanted_cell_anapri_data_indices_ub[
	missing_indices3
].astype(int)

full_wanted_cell_data_indices_lb = np.concatenate(
	(step2_wanted_cell_data_indices_lb, missing_anapri_cell_values_lb)
)
full_wanted_cell_data_indices_ub = np.concatenate(
	(step2_wanted_cell_data_indices_ub, missing_anapri_cell_values_ub)
)

full_wanted_cell_data_indices_lb = np.sort(full_wanted_cell_data_indices_lb)
full_wanted_cell_data_indices_ub = np.sort(full_wanted_cell_data_indices_ub)

full_wanted_cell_data_indices_lb = sorted(set(full_wanted_cell_data_indices_lb))
full_wanted_cell_data_indices_ub = sorted(set(full_wanted_cell_data_indices_ub))

# %%
# Create a dictionary of stylised cells
one_cell_index = []

stylised_cells_dict = {}

g = 0

for lb, ub in zip(full_wanted_cell_data_indices_lb, full_wanted_cell_data_indices_ub):
	lb = int(lb)
	ub = int(ub)

	m_ = full_data_indices[(full_data_indices > lb) & (full_data_indices < ub)].astype(
		int
	)
	stop_ = stop_indices[(stop_indices > lb) & (stop_indices < ub)].astype(int)
	g_neg = global_negator_indices[
		(global_negator_indices > lb) & (global_negator_indices < ub)
	].astype(int)
	pres_neg = precise_negator_indices[
		(precise_negator_indices > lb) & (precise_negator_indices < ub)
	].astype(int)
	pres_pos = precise_posator_indices[
		(precise_posator_indices > lb) & (precise_posator_indices < ub)
	].astype(int)
	b_ = inverse_term_indices[
		(inverse_term_indices > lb) & (inverse_term_indices < ub)
	].astype(int)
	mg_ = margin_term_indices[
		(margin_term_indices > lb) & (margin_term_indices < ub)
	].astype(int)
	a_ = full_anapri_data_indices[
		(full_anapri_data_indices > lb) & (full_anapri_data_indices < ub)
	].astype(int)
	map_ = mapping_term_indices[
		(mapping_term_indices > lb) & (mapping_term_indices < ub)
	].astype(int)
	or_ = ambivalence_term_indices[
		(ambivalence_term_indices > lb) & (ambivalence_term_indices < ub)
	].astype(int)
	b_ind = benign_indicator_term_indices[
		(benign_indicator_term_indices > lb) & (benign_indicator_term_indices < ub)
	].astype(int)
	pre_ind = premalignant_indicator_term_indices[
		(premalignant_indicator_term_indices > lb)
		& (premalignant_indicator_term_indices < ub)
	].astype(int)
	diff_ = differential_term_indices[
		(differential_term_indices > lb) & (differential_term_indices < ub)
	].astype(int)
	u_m = unspecified_malignancy_term_indices[
		(unspecified_malignancy_term_indices > lb)
		& (unspecified_malignancy_term_indices < ub)
	].astype(int)

	if m_.size > 0:
		for i in range(len(m_)):
			one_cell_index.append(m_[i])
	if a_.size > 0:
		for i in range(len(a_)):
			one_cell_index.append(a_[i])
	if u_m.size > 0:
		for i in range(len(u_m)):
			one_cell_index.append(u_m[i])
	if g_neg.size > 0:
		for i in range(len(g_neg)):
			one_cell_index.append(g_neg[i])

	cell_nucleus = np.sort(
		np.concatenate(
			(
				m_,
				stop_,
				g_neg,
				pres_neg,
				pres_pos,
				b_,
				mg_,
				a_,
				map_,
				or_,
				b_ind,
				pre_ind,
				diff_,
				u_m,
			)
		)
	)

	x = []

	for i in cell_nucleus:
		i = int(i)
		if list(cell_nucleus).count(i) > 1:
			if i in cat_data_indices_dict["Cs"]:
				if i in cat_data_indices_dict["Sa"]:
					if i not in a_:
						cell_nucleus = list(cell_nucleus)
						cell_nucleus.remove(i)
						cell_nucleus = np.array(cell_nucleus)

	for i in cell_nucleus:
		i = int(i)
		if i in m_:
			for key, value in cat_data_indices_dict.items():
				if np.any(np.in1d(value, i)):
					x.append(key)

		if i in a_:
			for key, value in roi_anapri_data_indices_dict.items():
				if np.any(np.in1d(value, i)):
					x.append(key)

		elif i in stop_:
			x.append(".")
		elif i in g_neg:
			x.append("Γ")
		elif i in pres_neg:
			x.append("¬")
		elif i in pres_pos:
			x.append("⌐")
		elif i in b_:
			x.append("B")
		elif i in mg_:
			x.append("∥")
		elif i in map_:
			x.append("⇔")
		elif i in or_:
			x.append("/")
		elif i in b_ind:
			x.append("β")
		elif i in pre_ind:
			x.append("δ")
		elif i in diff_:
			x.append("~")
		elif i in u_m:
			x.append("∪")

	x = [str(i) for i in x]

	stylised_cell = " ".join(x)

	stylised_cells_dict[g] = stylised_cell

	g = g + 1

# %%
# These are the coord_ticks closest to the terms
# wanted_coord_ticks = np.array(list(sorted(set(take_closest(coord_ticks, i) for i in one_cell_index))))
wanted_coord_ticks = np.array(
	list(sorted((take_closest(coord_ticks, i) for i in one_cell_index)))
)

# These are the indices of the coord_ticks, which will tell us the column
# the term came from
col_indices = np.where(np.in1d(coord_ticks, wanted_coord_ticks))[0]

# This tells us which column the malignant term reference came from
column_with_term = []
for w in col_indices:
	for key in column_map:
		if w in column_map[key]:
			column_with_term.append(key)
			break  # Break immediately we find the column

# This tells us which row the malignant term reference came from
row_ends = column_map[df.columns[df.shape[1] - 1]]
row_ticks = np.array(list(sorted(take_closest(row_ends, i) for i in col_indices)))

row_indices = []
for i in row_ticks:
	r = np.where(row_ends == i)[0][0]
	row_indices.append(r)

row_indices = list(np.array(row_indices) - 1)

# %%
stylised_cells = list(stylised_cells_dict.values())

interim_df = pd.DataFrame(
	{"row": row_indices, "column": column_with_term, "stylised cell": stylised_cells}
)

# Assuming you already have interim_df with the necessary columns
# Group by the 'row' column
grouped_df = interim_df.groupby("row")

# Initialize an empty dictionary to store the results
stylised_dict = {}

# Iterate over each group
for r, group in grouped_df:
	# Create a nested dictionary for each row
	stylised_dict[r] = {}

	# Iterate over the rows within the group
	for index, row in group.iterrows():
		stylised_dict[r][row["column"]] = row["stylised cell"]
# Now stylised_dict contains the desired mapping of column value to stylised cell value for each row

# %%
# with open(data_folder / "Processed/stylised_hp.txt", "w", encoding="utf-8") as f:
# 	for key, value in stylised_dict.items():
# 		f.write("%s: %s\n" % (int(key), value))
# f.close()

# %%
end = time.time()
print("The program took: ", (end - start) / 3600, "hrs", "to run.")
# 2.395 

# %%
