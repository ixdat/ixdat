from ixdat.techniques import ECMeasurement

path_to_file = "../test_data/biologic/" "Pt_poly_cv.mpt"

m = ECMeasurement.read(
    path_to_file,
    reader="biologic",
    name="ec_tools_test",
)

# ax = m.plot()

m.save()
i = m.id
del m

m1 = ECMeasurement.get(i)

m1.plot()
