from PyFoam.Infrastructure.CTestRun import CTestRun


class PlainIcoFoamCavity(CTestRun):
    def init(self):
        self.setParameters(solver="icoFoam",
                           originalCase="$FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity/",
                           sizeClass="tiny")


PlainIcoFoamCavity().run()
