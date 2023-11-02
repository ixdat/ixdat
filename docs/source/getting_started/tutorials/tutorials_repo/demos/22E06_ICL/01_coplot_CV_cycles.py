from ixdat import Measurement


ec = Measurement.read(
    "../../data/01/iridium_butterfly_short_CVA.mpt", reader="biologic"
)

ec.plot()

ec.calibrate(RE_vs_RHE=0.720, A_el=0.196)
cv = ec.as_cv()
cv.plot()

cv.plot_cycles()

ax = cv[1].plot(color="k", label="early")
cv[200].plot(color="b", ax=ax, label="late")
ax.legend()

cv.redefine_cycle(start_potential=0.5, redox=True)

ax = cv[1].plot(color="k", label="early")
cv[200].plot(color="b", ax=ax, label="late")
ax.legend()
