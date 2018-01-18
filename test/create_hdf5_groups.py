import h5py as h5

# matlab hdf5 matrix, CSC
f = h5.File("matrix_00.mat", "r")
fexp = f["Aexp"]
ir = fexp["ir"]
jc = fexp["jc"]
aij = fexp["data"]

# scipy sparse, CSC
# nrow x nrow
# nnz
# nexp
mat_csc = scipy.sparse.csc_matrix((aij, (ir, jc)), shape=(nrow, nrow), dtype=">f8", copy=False)

Ir_pre
profilesFile=fullfile(matrixPath,'Data','profile_data');

# scipy sparse, CSR
mat_csr = mat.tocsr()

# for storage, petsc standard view format
nnzrow = numpy.diff(mat_csr.indptr)
colidx = mat_csr.indices
aij = mat_csr.data

f = h5.File("test_group.nc", "w")
dset = f.create_dataset("nnzrow", (nrow,), dtype = ">i4", data=nnzrow)#, fillvalue=None, compression="gzip")
dset = f.create_dataset("colidx", (nnz,), dtype = ">i4", data=colidx)#, fillvalue=None, compression="gzip")
dset = f.create_dataset("aij", (nnz,), dtype = ">f8", data=aij)#, fillvalue = None, compression = "gzip")


# for creation, petsc mataij
nv = nrow, assert
nvloc
nexp
nnz
indptr = numpy.append(0, numpy.cumsum(nnzrow))
colidx = colidx
for ie in nexp:
    A = PETSc.Mat()
    A.setType('aij') # PETSc.Mat.Type.AIJ
    A.setSize(((nvloc, nv), (nvloc, nv)))
    A.setPreallocationNNZ(nnz)
    A.setValuesCSR(indptr, colidx, aij[ie])
    A.assemble()




for i in range(nrow):
    A[i, colidx[indptr[i]:indptr[i+1]]] = aij[indptr[i]:indptr[i+1]]
A.setValues(self, rows, cols, values, addv=None)

nnzcol = numpy.diff(jc)
A.setPreallocationNNZ(nnz)
for j in range(N):
    A[ir[jc[j]:jc[j+1]], j] = aij[jc[j]:jc[j+1]]

for i in range(N):
    A[i, jc[ir[i]:ir[i+1]]] = aij[ir[i]:ir[i+1]]

#group = f.create_group("Aexp_0")
#dset = group.create_dataset("ir", (10,), dtype = ">i4", fillvalue = None, compression = "gzip")
#dset = group.create_dataset("jc", (100,), dtype = ">i4", fillvalue = None, compression = "gzip")
#dset = group.create_dataset("aij", (100,), dtype = ">f8", fillvalue = None, compression = "gzip")

#group = f.create_group("Aexp_1")
#dset = group.create_dataset("ir", (10,), dtype = ">i4", fillvalue = None, compression = "gzip")
#dset = group.create_dataset("jc", (100,), dtype = ">i4", fillvalue = None, compression = "gzip")
#dset = group.create_dataset("aij", (100,), dtype = ">f8", fillvalue = None, compression = "gzip")

#dispatch_key(f, conf_key, conf_key_dict)
#create_attributes(f, conf_dict["Global attribute"])

f.close()

