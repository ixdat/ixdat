from ixdat.techniques import ECMeasurement

path_to_file = (
    "../test_data/biologic_mpt_and_zilien_tsv/"
    "2020-07-29 10_30_39 Pt_poly_cv_01_02_CVA_C01.mpt"
)

m = ECMeasurement.read(path_to_file, reader="biologic", name="ec_tools_test",)

# ax = m.plot()

m.save()
i = m.id
del m

m1 = ECMeasurement.get(i)

m1.plot()
