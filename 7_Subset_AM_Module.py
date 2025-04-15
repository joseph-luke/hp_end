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
file_to_open_1 = data_folder / "Processed/am_basis.csv"
df_basis = pd.read_csv(file_to_open_1, low_memory=False)
basis = df_basis.iloc[0, 0]
B_synthetic = [basis.split()]

# %%
data = {}
with open(data_folder / "Processed/condensed_am.txt", "r", encoding="utf-8") as file:
	for line in file:
		s = line.split(":", maxsplit=1)
		x = int(s[0])
		t = s[1].split("\n")[0].strip()
		data[x] = t
file.close()

# %%
space = list(data.values())
T_test = [t.split() for t in space]

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
symbols.append("m")
symbols.append(".")
symbols.append("|")

codes = [str(int) for int in sieve[1: (len(symbols) + 1)]]

nthumba_dict = dict(zip(symbols, codes))
inv_nthumba_dict = {v: k for k, v in nthumba_dict.items()}

# %%
T_hat_test = []
for t_ in T_test:
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
file_to_save = data_folder / "Processed/predict_am.csv"
fit_df.to_csv(file_to_save, index=False)

# %%
# Finds best basis + best basis model for next step
import rpy2.robjects as robjects

robjects.r.source("subset_am.R")

# %%
file_to_open_1 = data_folder / "Processed/condensed_am.csv"
df = pd.read_csv(file_to_open_1, low_memory=False)

file_to_open_1 = data_folder / "Processed/potential_am.csv"
df_potential = pd.read_csv(file_to_open_1, low_memory=False)

file_to_open_1 = data_folder / "Processed/predicted_am.csv"
df_predictions = pd.read_csv(file_to_open_1, low_memory=False)

df_potential['prediction'] = df_predictions.iloc[:,0]

# %%
grouped_df = df_potential.groupby('record', sort=False)

# %%
# Initialize an empty dictionary to store the results
mapped_dict = {}
cleaned_mapped_dict = {}

# Iterate over each group
for r, group in grouped_df:
	row_primes = []
	row_primes_to_delete = []

	no_row_primes = []

	for index, row in group.iterrows():
		if row['prediction'] == 1:
			row_primes.append(int(row['anatomyprime']))
		if sum(group.prediction) == 0:
			no_row_primes.append(row['malignancy'])

	row_primes = sorted(set(row_primes))
	no_row_primes = sorted(set(no_row_primes))

	for inx, val in enumerate(row_primes):
		val = int(val)
		if inx < len(row_primes)-1:
			for i in range(inx+1, len(row_primes)):
				if val % row_primes[i] == 0:
					row_primes_to_delete.append(row_primes[i])
				elif row_primes[i] % val == 0:
					row_primes_to_delete.append(val)

	# Create a nested dictionary for each row
	mapped_dict[r] = {}
	cleaned_mapped_dict[r] = ""

	for index, row in group.iterrows():
		if row['prediction'] == 1:
			if row['anatomyprime'] not in row_primes_to_delete:
				mapped_dict[r][row['malignancy']] = []

	for index, row in group.iterrows():
		if row['prediction'] == 1:
			if row['anatomyprime'] not in row_primes_to_delete:
				mapped_dict[r][row['malignancy']].append(f"{row['anatomy']}".strip())

	items = []
	for key, value in mapped_dict[r].items():
		if value:
			items.append(f" {key} ⬌ {", ".join(value)} ")

	cleaned_mapped_dict[r] += "; ".join(items)
	
	if no_row_primes:
		cleaned_mapped_dict[r] = f"{"; ".join(sorted(set(no_row_primes)))} inferred"
		
	cleaned_mapped_dict[r] = cleaned_mapped_dict[r].strip()

# %%
condensed_map_dict = {}

for key, value in cleaned_mapped_dict.items():
	new_value = value
	if new_value:
		condensed_map_dict[key] = new_value
	if not new_value:
		condensed_map_dict[key] = 'No specific map found'

# %%
df_am = df.iloc[list(mapped_dict.keys())]

df_am.loc[:, "A-M Map"] = list(condensed_map_dict.values())

file_to_save = data_folder / "Processed/am_records.xlsx" # For clean viewing
df_am.to_excel(file_to_save, index=False, engine='xlsxwriter')

# %%
