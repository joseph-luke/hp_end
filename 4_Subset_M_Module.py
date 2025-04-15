# %%
import numpy as np
from pathlib import Path
import pandas as pd
import re
from sympy import sieve

data_folder = Path(
	"C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"
)
# %%
# Import best basis
file_to_open_1 = data_folder / "Processed/m_basis.csv"
df_basis = pd.read_csv(file_to_open_1, low_memory=False)
basis = df_basis.iloc[0, 0]
B_synthetic = [basis.split()]

# %%
data = {}
with open(data_folder / "Processed/condensed_m.txt", "r", encoding="utf-8") as file:
	for line in file:
		s = line.split(":", maxsplit=1)
		x = int(s[0])
		t = s[1].split("\n")[0].strip()
		data[x] = t
file.close()

file_to_open_1 = data_folder / "Processed/condensed_m.csv"
df = pd.read_csv(file_to_open_1)

# %%
test_list = list(data.keys())
pre_space = [data[i] for i in test_list]
space = []

for vec in pre_space:
	vec = re.split(r"\.", vec)
	for v in vec:
		v = v.strip()
		symbol_list = list(re.findall(r"\w+|[^\w\s]", v))
		space.append(symbol_list)

# %%
v_record = []

n = 0
for inx, vec in zip(test_list, pre_space):
	vec = re.split(r"\.", vec)
	vec2 = []
	for v in vec:
		vec2.append(v.strip())
	for v in vec2:
		v_record.append(inx)

# %%
symbols = []
symbols.append("Γ")
symbols.append("¬")
symbols.append("⌐")
symbols.append("B")
symbols.append("∥")
symbols.append("⇔")
symbols.append("/")
symbols.append("β")
symbols.append("δ")
symbols.append("~")
symbols.append("A")
symbols.append("M")

codes = [str(int) for int in sieve[1 : (len(symbols) + 1)]]

nthumba_dict = dict(zip(symbols, codes))
inv_nthumba_dict = {v: k for k, v in nthumba_dict.items()}

# %%
T_hat_test = []
for t_ in space:
	t_hat_test = list((pd.Series(t_)).map(nthumba_dict))
	T_hat_test.append(t_hat_test)

T_test = T_hat_test

# %%
# Fit b_ to T_
operations = []

for basis in B_synthetic:
	b_pristine = list((pd.Series(basis)).map(nthumba_dict))

	T_ = [t_ for t_ in T_test]

	for t_ in T_:
		t_pristine = t_.copy()

		b_ = np.array(b_pristine)

		t_ = np.array(t_pristine)

		# b_ space
		# Tabulate operations
		b_t_operations = []
		b_t_matrix = []
		for tt in t_:
			col = []
			for bb in b_:
				if tt == bb:
					col.append(1)
				elif tt != bb:
					col.append(0)
			b_t_matrix.append(col)

		b_t_matrix = np.array(b_t_matrix).T

		diagonal_y = list(np.arange(b_t_matrix.shape[0]))
		diagonal_x = list(np.arange(b_t_matrix.shape[1]))
		diagonal_yx = np.intersect1d(diagonal_y, diagonal_x)
		diagonal_yx = np.vstack((diagonal_yx, diagonal_yx)).T

		search_space = np.indices(b_t_matrix.shape).reshape(2, -1).T
		ones_indices = np.array(np.where(b_t_matrix == 1)).reshape(2, -1).T
		for diagonal in diagonal_yx:
			# yx is the value on the diagonal
			yx = b_t_matrix[diagonal[0], diagonal[1]]

			if yx == 1:
				# Delete the diagonal column and row from the search space
				search_space = np.delete(
					search_space, np.where((search_space == diagonal[0]).any(1))[0], 0
				)
				# Delete the diagonal column and row from the ones indices
				ones_indices = np.delete(
					ones_indices, np.where((ones_indices == diagonal).any(1))[0][0], 0
				)
				pass

			elif yx != 1:
				available_ones = np.array(
					[
						x
						for x in set(tuple(x) for x in search_space)
						& set(tuple(x) for x in ones_indices)
					]
				)

				# Sorts based on row
				if available_ones.shape[0] > 1:
					available_ones = available_ones[np.argsort(available_ones[0, :])]

				if available_ones.shape[0] > 0:
					nearest_one = available_ones[0]
					row_distance = nearest_one[0] - diagonal[0]
					col_distance = nearest_one[1] - diagonal[1]

					row_flips = np.flip(np.arange(row_distance + 1))
					col_flips = np.flip(np.arange(col_distance + 1))

					for r in row_flips:
						if r >= 1:
							b_t_matrix[[diagonal[0] + r, diagonal[0] + r - 1], :] = (
								b_t_matrix[[diagonal[0] + r - 1, diagonal[0] + r], :]
							)
							# Rearrange b_ as well
							b_indx = np.arange(len(b_))
							b_indx2 = b_indx.copy()
							b_indx2[diagonal[0] + r - 1] = b_indx[diagonal[0] + r]
							b_indx2[diagonal[0] + r] = b_indx[diagonal[0] + r - 1]
							b_ = b_[b_indx2]
							b_t_operations.append(
								int(b_[diagonal[0] + r - 1]) * int(b_[diagonal[0] + r])
							)

					for c in col_flips:
						if c >= 1:
							b_t_matrix[:, [diagonal[1] + c, diagonal[1] + c - 1]] = (
								b_t_matrix[:, [diagonal[1] + c - 1, diagonal[1] + c]]
							)
							# Rearrange t_ as well
							t_indx = np.arange(len(t_))
							t_indx2 = t_indx.copy()
							t_indx2[diagonal[1] + c - 1] = t_indx[diagonal[1] + c]
							t_indx2[diagonal[1] + c] = t_indx[diagonal[1] + c - 1]
							t_ = t_[t_indx2]
							b_t_operations.append(
								int(t_[diagonal[1] + c - 1]) * int(t_[diagonal[1] + c])
							)

					# Re-index ones
					ones_indices = np.array(np.where(b_t_matrix == 1)).reshape(2, -1).T

					# Delete the diagonal column and row from the search space
					search_space = np.delete(
						search_space,
						np.where((search_space == diagonal[0]).any(1))[0],
						0,
					)
					# Delete the diagonal column and row from the ones indices
					ones_indices = np.delete(
						ones_indices,
						np.where((ones_indices == diagonal).all(1))[0][0],
						0,
					)

					available_ones = np.array(
						[
							x
							for x in set(tuple(x) for x in search_space)
							& set(tuple(x) for x in ones_indices)
						]
					)

		# Shifts
		# We shift on whichever axis is longer
		n = b_t_matrix.shape[0]
		m = b_t_matrix.shape[1]

		row_shift = np.arange(n)
		col_shift = np.arange(m)

		for diagonal in diagonal_yx:
			# yx is the value on the diagonal
			yx = b_t_matrix[diagonal[0], diagonal[1]]
			if yx == 1:
				row_shift = np.delete(row_shift, np.where(row_shift == diagonal[0]))
				col_shift = np.delete(col_shift, np.where(col_shift == diagonal[1]))

			elif yx != 1:
				# Retain the row and column containing the non-1 diagonal
				pass

		b_t_matrix = np.delete(b_t_matrix, row_shift, 0)
		b_t_matrix = np.delete(b_t_matrix, col_shift, 1)

		for shift in row_shift:
			b_t_operations.append(-int(b_[shift]))

		for shift in col_shift:
			b_t_operations.append(-int(t_[shift]))

		sq_reduced_operations = []
		for operation in b_t_operations:
			if not np.sqrt(abs(operation)) % 1 == 0:
				sq_reduced_operations.append(operation)

		operations.append(list(sorted(set(b_t_operations), reverse=True)))

# %%
operators_and_shifts = np.unique(np.concatenate(operations))[::-1]

fit = np.vstack(
	[np.isin(operators_and_shifts, operation) for operation in operations]
).astype(int)
fit_df = pd.DataFrame(fit)

fit_df.columns = operators_and_shifts

# %%
file_to_save = data_folder / "Processed/predict_m.csv"
fit_df.to_csv(file_to_save, index=False)

# %%
# Finds best basis + best basis model for next step
import rpy2.robjects as robjects

robjects.r.source("subset_m.R")

# %%
# Import basis predictions
file_to_open_1 = data_folder / "Processed/predicted_m.csv"
df_predictions = pd.read_csv(file_to_open_1, low_memory=False)

# %%
set_v_records = list(sorted(set(v_record)))


m_predictions_df = pd.DataFrame(
	{"record": v_record, "prediction": df_predictions.iloc[:, 0]}
)

# Group by the 'row' column
grouped_df = m_predictions_df.groupby("record")

# Initialize an empty dictionary to store the results
record_M_state = {}

# Iterate over each group
for r, group in grouped_df:
	# Create a nested dictionary for each row
	record_M_state[r] = {}

	# Iterate over the rows within the group
	for index, row in group.iterrows():
		record_M_state[r] = row["prediction"]

total_m_score = []

for v in record_M_state.values():
	if v < 1:
		total_m_score.append("Hypothesised Other Record")
	elif v > 0:
		total_m_score.append("Hypothesised Malignant Record")


# %%
df_all_predicitions = df.iloc[set_v_records]
df_all_predicitions.loc[:, "SpecimenClass"] = total_m_score

# For moving on to the next timeline steps:
df_m = df_all_predicitions.loc[
	df_all_predicitions.SpecimenClass == "Hypothesised Malignant Record"
]

file_to_save = data_folder / "Processed/m_records.csv"
df_m.to_csv(file_to_save, index=False)

# %%
