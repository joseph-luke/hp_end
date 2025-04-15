# %%
import itertools
import numpy as np
from pathlib import Path
import pandas as pd
from sympy import sieve


data_folder = Path(
	"C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"
)

# %%
M_ = [
	["A", "M", "⌐"],
	["A", "M"],
	["A", "⇔", "A", "M"],
	["M"],
	["⇔", "M"],
	["⇔", "M", "A"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M", "∥"],
	["M", "A"],
	["A", "⇔", "M", "/"],
	["A", "⇔", "A", "M", "¬", "A"],
	["A", "⇔", "M", "A", "∥"],
	["A", "M", "δ"],
	["A", "β", "⇔", "M"],
	["⇔", "A", "M"],
	["A", "⇔", "M", "/", "A"],
	["A", "⇔", "M", "~", "A"],
	["M", "~", "M"],
	["A", "¬", "⇔", "M"],
	["⇔", "A", "M", "⇔", "A", "∥", "⇔", "A"],
	["A", "M", "A", "M"],
	["M", "A", "~", "M"],
	["A", "⇔", "β", "M"],
	["M", "¬"],
	["M", "A", "M"],
	["A", "β", "⇔", "A", "M", "⇔"],
	["M", "~", "A", "M"],
	["A", "⇔", "M", "A", "δ"],
	["A", "⇔", "M", "∥", "M", "δ"],
	["A", "⇔", "M", "⇔", "δ"],
	["A", "⇔", "M", "~", "M"],
	["A", "⇔", "M", "Γ"],
	["A", "⇔", "M", "A", "/", "A"],
	["A", "⇔", "M", "¬"],
	["A", "⇔", "M", "/", "~"],
	["A", "⇔", "M", "¬", "∥"],
	["A", "⇔", "M", "~", "A", "/", "A"],
	["A", "⇔", "A", "M", "⇔", "A", "M"],
	["A", "⇔", "β", "⇔", "M"],
	["A", "M", "A", "/", "A"],
	["A", "¬", "⇔", "A", "⇔", "M", "A"],
	["β", "⇔", "A", "M"],
	["β", "M"],
	["⇔", "β", "M"],
	["⇔", "M", "∥", "¬"],
	["M", "⇔"],
	["M", "⌐"],
	["M", "⇔", "M"],
]

O_ = [
	["¬", "A", "M", "⌐"],
	["A", "⇔", "β", "M", "B"],
	["A", "⇔", "¬", "A", "M", "⇔"],
	["¬", "M"],
	["A", "⇔", "B", "β", "M"],
	["A", "M", "⇔", "B", "¬"],
	["⇔", "¬", "M"],
	["⇔", "Γ", "M"],
	["A", "⇔", "M", "¬", "⌐"],
	["A", "⌐", "Γ", "A", "M"],
	["⇔", "¬", "M", "/"],
	["A", "β", "⇔", "B", "M"],
	["A", "⇔", "B", "M"],
	["A", "⇔", "β", "M"],
	["A", "⇔", "A", "⇔", "M", "δ"],
	["β", "M", "β", "A"],
	["⇔", "A", "¬"],
	["β", "M", "¬"],
	["M", "δ", "A"],
	["Γ", "M"],
	["M", "β", "M", "A"],
	["A", "⇔", "¬", "M"],
	["B", "A", "δ", "Γ", "M", "A"],
	["A", "⇔", "M", "δ", "⇔", "¬", "⇔", "∥"],
	["A", "⇔", "M", "δ"],
	["A", "⇔", "δ", "M"],
	["A", "β", "M", "B"],
	["A", "β", "⇔", "B", "β", "M"],
	["⇔", "¬", "M", "Γ"],
	["β", "M", "B"],
	["B", "A", "Γ", "M"],
	["δ", "M"],
	["δ", "/", "M"],
	["δ", "M", "⇔", "A"],
	["Γ", "A", "/", "M"],
	["M", "δ"],
	["M", "Γ"],
]

M_symbolic = M_.copy()
O_symbolic = O_.copy()

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
M_hat = []
for m_ in M_:
	if m_ not in O_:
		m_hat = list((pd.Series(m_)).map(nthumba_dict))
		M_hat.append(m_hat)

O_hat = []
for o_ in O_:
	if o_ not in M_:
		o_hat = list((pd.Series(o_)).map(nthumba_dict))
		O_hat.append(o_hat)

M_ = M_hat
O_ = O_hat

B_ = [b_ for b_ in M_]
B_ += [b_ for b_ in O_]

# %%
# # Fit b_ to T_
# b_symbols = ['A', '⇔', 'δ', 'M', 'Γ', 'β', '¬']
# B_synthetic = list(itertools.permutations(b_symbols))

B_synthetic = [["M", "δ", "Γ", "β", "⇔", "A", "¬"]]
B_synthetic += [["M", "Γ", "β", "⇔", "A", "¬"]]

B_synthetic = [list(i) for i in B_synthetic]
B_synthetic_symbolic = B_synthetic.copy()

# %%
operations = []

basis_used = []

y = []

for basis in B_synthetic:
	b_pristine = list((pd.Series(basis)).map(nthumba_dict))

	T_ = [t_ for t_ in B_ if t_ != b_pristine]

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

		if list(t_pristine) in M_:
			basis_used.append(f"{" ".join([inv_nthumba_dict[b] for b in b_pristine])}")
			operations.append(list(sorted(set(b_t_operations), reverse=True)))
			y.append(1)

		elif list(t_pristine) in O_:
			basis_used.append(f"{" ".join([inv_nthumba_dict[b] for b in b_pristine])}")
			operations.append(list(sorted(set(b_t_operations), reverse=True)))
			y.append(0)

# %%
operators_and_shifts = np.unique(np.concatenate(operations))[::-1]

fit = np.vstack(
	[np.isin(operators_and_shifts, operation) for operation in operations]
).astype(int)
fit_df = pd.DataFrame(fit)

fit_df.columns = operators_and_shifts
fit_df["y"] = y
fit_df.insert(0, "basis", basis_used)

# %%
B_synthetic = B_synthetic_symbolic.copy()
B_synthetic_hat = [list((pd.Series(b_)).map(nthumba_dict)) for b_ in B_synthetic]
B_synthetic = B_synthetic_hat
B_ = [b_ for b_ in B_synthetic]

# %%
# Test b_ performance
O_test = [
	["B", "A", "¬", "A", "M", "/", "M"],
	["¬", "M"],
	["A", "⇔", "β", "M", "β"],
	["A", "⇔", "A", "¬", "M"],
	["⇔", "M", "δ"],
	["¬", "A", "δ", "/", "A", "M"],
	["β", "M", "β"],
	["Γ", "/", "M", "⇔", "A"],
	["⇔", "¬", "M"],
	["¬", "/", "M"],
	["Γ", "A", "M"],
	["A", "⇔", "A", "¬", "M"],
	["¬", "M", "/"],
	["A", "⇔", "¬", "M"],
	["Γ", "A", "M"],
	["¬", "A", "M"],
	["¬", "/", "M"],
	["¬", "M"],
	["~", "M"],
	["~", "A", "β", "A", "B", "A", "M"],
	["⌐", "A", "~", "A", "M"],
	["M", "⌐", "⇔"],
	["¬", "M", "/"],
	["Γ", "A", "M"],
	["Γ", "M"],
	["¬", "A", "M", "/", "M"],
	["M", "δ"],
	["⇔", "M"],
	["A", "⇔", "¬", "M"],
	["A", "⇔", "~", "M"],
	["A", "⇔", "⌐", "¬", "/", "M"],
	["⇔", "¬", "M"],
	["M", "δ"],
	["∥", "¬", "δ", "M"],
	["∥", "⇔", "∥", "δ", "M"],
	["¬", "~", "M"],
	["A", "⇔", "¬", "M"],
	["⌐", "¬", "δ", "/", "M"],
	["A", "⇔", "¬", "M"],
	["A", "⇔", "¬", "/", "M"],
	["A", "⇔", "M", "δ"],
	["¬", "M", "∥", "¬", "/", "A", "¬", "M"],
	["M", "⌐", "~", "A", "¬", "~"],
	["A", "Γ", "δ", "M", "/", "A", "M"],
	["¬", "M", "/"],
	["¬", "M", "/"],
	["¬", "M"],
	["A", "⇔", "β", "M"],
	["¬", "M"],
	["¬", "/", "M"],
	["¬", "M"],
	["A", "⇔", "B", "β", "M"],
	["A", "⇔", "M", "B"],
	["A", "⌐", "Γ", "M"],
	["A", "⇔", "M", "δ"],
	["A", "⇔", "M", "Γ"],
	["β", "M", "⇔", "A"],
	["¬", "/", "M"],
	["Γ", "A", "M", "/", "M"],
	["M", "¬"],
	["¬", "M"],
	["¬", "M"],
	["δ", "M", "A"],
	["⇔", "Γ", "M"],
	["A", "⇔", "B", "β", "A", "¬", "δ", "/", "M"],
	["A", "¬", "⇔", "∥", "M"],
	["¬", "M", "/"],
	["Γ", "A", "M"],
	["M", "δ"],
	["⌐", "¬", "M"],
	["¬", "M"],
	["¬", "A", "M", "¬", "⌐"],
	["~", "M", "⇔", "⌐", "¬", "⇔"],
	["A", "⇔", "A", "¬", "M"],
	["⇔", "~", "M", "⇔", "A", "/", "M", "A", "/", "A"],
	["A", "⇔", "A", "Γ", "/", "M"],
	["¬", "M"],
	["A", "⇔", "M", "δ", "⇔", "¬", "M", "⇔", "A", "∥", "¬"],
	["¬", "/", "A", "M", "A"],
	["A", "⇔", "¬", "M"],
	["A", "⇔", "M", "δ", "¬"],
	["M", "δ"],
]

M_test = [
	["A", "⇔", "M", "A", "∥"],
	["A", "⇔", "M"],
	["M", "A"],
	["⇔", "⌐", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["β", "M", "⌐", "~", "M", "/", "~", "A", "/", "A"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["⌐", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M", "⇔"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["⇔", "M", "⌐", "A", "⇔"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["⇔", "M"],
	["M", "A", "M", "A"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "A", "M", "⇔", "A", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "M"],
	[
		"⇔",
		"A",
		"⇔",
		"¬",
		"⇔",
		"M",
		"⇔",
		"A",
		"⇔",
		"A",
		"⇔",
		"∥",
		"⇔",
		"∥",
		"A",
		"⇔",
		"A",
		"⇔",
		"A",
	],
	["M"],
	["A", "⇔", "M", "⇔", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "⇔", "A"],
	["M"],
	["M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "M", "⇔", "A"],
	["A", "⇔", "M", "⇔"],
	["A", "⇔", "A", "M"],
	["⇔", "M"],
	["⇔", "M", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M", "/", "M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["M"],
	["⇔", "A", "⇔", "M", "⇔", "A", "⇔", "∥", "⇔", "A", "⇔", "A", "⇔"],
	["M", "A"],
	["A", "⇔", "M", "⇔", "A"],
	["A", "⇔", "M"],
	["M", "⇔", "A"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "A"],
	["⇔", "A", "M", "~"],
	["M"],
	["A", "⇔", "M", "A"],
	["M"],
	["A", "⇔", "M"],
	["M", "A"],
	["M", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["M", "A"],
	["⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "⇔", "M", "⇔", "A"],
	["⇔", "M", "A", "⌐", "A"],
	["⇔", "A", "M"],
	["⇔", "M", "/"],
	["M"],
	["A", "⇔", "A", "M"],
	["⇔", "M"],
	["A", "⇔", "A", "M", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "⇔", "M"],
	["A", "⇔", "M", "A", "∥"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["M"],
	["A", "⇔", "M"],
	["A", "M"],
	["M", "A"],
	["A", "⇔", "M", "A"],
	["⇔", "A", "M", "⇔", "A", "⇔", "¬", "A", "⇔", "A", "M"],
	["A", "⇔", "M", "A", "⇔", "A", "⇔", "∥", "⇔", "∥"],
	["A", "⇔", "A", "M"],
	["⇔", "M"],
	["A", "⇔", "M", "A", "M"],
	["A", "M", "/", "A", "M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "M", "⇔", "A", "M"],
	["⇔", "M", "⇔", "A"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["δ", "M", "⇔", "A", "/", "M"],
	["A", "⇔", "M", "A"],
	["⇔", "A", "M"],
	["A", "⇔", "A", "M"],
	["M", "A"],
	["⇔", "A", "M", "A", "M"],
	["M", "⇔", "A"],
	["M", "A"],
	["A", "M"],
	["M", "A", "⌐"],
	["M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "A", "M"],
	["⇔", "A", "M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "Γ", "M"],
	["A", "⇔", "A", "M", "⇔"],
	["⇔", "M"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "∥", "A"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "M"],
	["A", "⇔", "M", "A", "⌐", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["⇔", "M"],
	["A", "⇔", "M", "¬", "A"],
	["A", "⇔", "M", "A"],
	["M"],
	["A", "⇔", "M"],
	["A", "M"],
	["⇔", "M"],
	["∥", "M", "⇔", "A", "∥"],
	["A", "⇔", "A", "M", "⇔"],
	["A", "⇔", "M", "A"],
	["A", "M"],
	["⇔", "A", "M", "B"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M", "¬", "A"],
	["A", "⇔", "M", "A", "/", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "M", "¬", "⇔", "A", "⇔", "A", "⇔", "A", "M", "⇔"],
	["A", "M"],
	["A", "⇔", "M", "A", "⇔", "A"],
	["M", "⇔", "A"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["M", "A"],
	["M"],
	["A", "⇔", "M", "⇔", "A"],
	["M", "A"],
	["A", "⇔", "M"],
	["M"],
	["A", "¬", "⇔", "M", "⇔", "¬"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "¬", "A"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "A"],
	["M", "A"],
	["¬", "M"],
	["⇔", "M"],
	["A", "M"],
	["M"],
	["A", "⇔", "M"],
	["M", "⇔", "A"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M", "A", "⇔"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "⇔", "M", "⇔", "A"],
	["A", "⇔", "M"],
	["M", "A"],
	["A", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M", "A", "M"],
	["M"],
	["A", "⇔", "M"],
	["M", "A"],
	["A", "⇔", "A", "M"],
	["M"],
	["A", "⇔", "M", "⇔", "A"],
	["A", "⇔", "A", "M"],
	["⇔", "M", "⇔", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["M"],
	["⇔", "M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["⇔", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "⇔", "A"],
	["⇔", "A", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "A", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M", "A", "M", "A", "M"],
	["A", "⇔", "M"],
	["A", "M", "A"],
	["M", "A"],
	["A", "⇔", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "⇔", "∥", "¬"],
	["M", "A"],
	["M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["A", "M", "A", "M"],
	["M"],
	["A", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "A", "⇔", "A", "M", "/", "A"],
	["A", "M"],
	["~", "A", "M", "A"],
	["M", "⇔", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["M", "A"],
	["A", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["M"],
	["A", "⇔", "M"],
	["⌐", "M"],
	["A", "M"],
	["A", "M", "A"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "δ", "M"],
	["⇔", "M", "⇔", "∥"],
	["M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["A", "⇔", "M"],
	["M", "A"],
	["M"],
	["M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "A", "/", "A", "/"],
	["M"],
	["A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["⇔", "M", "A"],
	["A", "⇔", "M"],
	["A", "M"],
	["⇔", "A", "M", "⇔", "A", "⇔", "¬", "A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "A", "M"],
	["A", "⇔", "M"],
	["M", "⇔", "A"],
	["A", "⇔", "M", "A"],
	["A", "M"],
	["M", "A"],
	["⇔", "A", "M", "A"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["M", "A", "⌐", "¬"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "δ", "M", "⇔", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["M", "β"],
	["~", "A", "⇔", "A", "M", "A", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "β", "M", "∥", "⌐"],
	["A", "⇔", "M"],
	["A", "⇔", "M", "A"],
	["A", "⇔", "M", "¬", "A"],
	["A", "⇔", "A", "M"],
	["M"],
	["M", "A"],
	["M", "¬", "A", "⌐", "A"],
	["M", "A"],
	["A", "⇔", "A", "M"],
	["M", "⇔", "A"],
	["M", "A"],
	["A", "⇔", "A", "δ", "M"],
	["M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M"],
	["M"],
	["M", "A"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["A", "⇔", "M"],
	["⇔", "M"],
	["M"],
	["A", "¬", "A", "⇔", "M", "A"],
	["A", "⇔", "M"],
	["A", "M"],
	["A", "⇔", "M", "⇔"],
	["A", "⇔", "M", "A"],
	["β", "M", "⇔", "A"],
]

# %%
M_hat_test = []
for m_ in M_test:
	if m_ not in O_test:
		m_hat_test = list((pd.Series(m_)).map(nthumba_dict))
		M_hat_test.append(m_hat_test)

O_hat_test = []
for o_ in O_test:
	if o_ not in M_test:
		o_hat_test = list((pd.Series(o_)).map(nthumba_dict))
		O_hat_test.append(o_hat_test)

M_test = M_hat_test
O_test = O_hat_test

T_test = M_test + O_test

# %%
# Fit data and record basis error
operations = []

basis_used = []

y = []

for b_ in B_:
	b_pristine = b_.copy()

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

		if list(t_pristine) in M_test:
			basis_used.append(f"{" ".join([inv_nthumba_dict[b] for b in b_pristine])}")
			operations.append(list(sorted(set(b_t_operations), reverse=True)))
			y.append(1)

		elif list(t_pristine) in O_test:
			basis_used.append(f"{" ".join([inv_nthumba_dict[b] for b in b_pristine])}")
			operations.append(list(sorted(set(b_t_operations), reverse=True)))
			y.append(0)

# %%
operators_and_shifts = np.unique(np.concatenate(operations))[::-1]

fit2 = np.vstack(
	[np.isin(operators_and_shifts, operation) for operation in operations]
).astype(int)
fit_df2 = pd.DataFrame(fit2)

fit_df2.columns = operators_and_shifts
fit_df2["y"] = y
fit_df2.insert(0, "basis", basis_used)

# %%
file_to_save = data_folder / "Processed/fit_m.csv"
fit_df.to_csv(file_to_save, index=False)

file_to_save = data_folder / "Processed/test_m.csv"
fit_df2.to_csv(file_to_save, index=False)

# %%
# Finds best basis + best basis model for next step
import rpy2.robjects as robjects
robjects.r.source("fit_m.R")

# %%
