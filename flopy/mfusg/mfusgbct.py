import numpy as np
from typing import List, Dict, Tuple, Union
from ..pakbase import Package
from ..utils import Util3d, Util2d
from ..utils.flopy_io import line_parse
from ..utils.utils_def import get_open_file_object
from ..utils.utils_def import get_util2d_shape_for_layer as get_u2d_shape
from .mfusg import MfUsg


class MfUsgBct(Package):
    """
    Block centered transport package class for MODFLOW-USG
    """

    def __init__(
        self,
        model: MfUsg,
        itrnsp: int = 1,
        ibctcb: int = 0,
        mcomp: int = 1,
        ic_ibound_flg: int = 1,
        itvd: int = 1,
        iadsorb: int = 0,
        ict: int = 0,
        cinact: float = -999.0,
        ciclose: float = 1.0e-6,
        idisp: int = 1,
        ixdisp: int = 0,
        diffnc: float = 0.0,
        izod: int = 0,
        ifod: int = 0,
        icbund: int = 1,
        ifmbc: int = 0,
        iheat: int = 0,
        imcomp: int = 0,
        idispcln: int = 0,
        nseqitr: int = 0,
        options: Union[List[str], None] = None,
        porosity: float = 0.1,
        bulkd: float = 2.6,
        dispd: Union[Dict[str, float], None] = None,
        adsorb: Union[List[float], None] = None,
        sconc: List[float] = 0.0,
        awadsorbd: Union[Dict[str, float], None] = None,
        extension: str = "bct",
        filenames: Union[str, None] = None,
        unitnumber: Union[int, None] = None,
        unitnumber_flowtransport: Union[Tuple[int], None] = None,
    ) -> None:
        # set default unit number of one is not specified
        if unitnumber is None:
            unitnumber = MfUsgBct._defaultunit()

        # set filenames
        filenames = self._prepare_filenames(filenames)

        # call base package constructor
        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        dis = model.get_package("DIS")
        if dis is None:
            dis = model.get_package("DISU")

        structured = model.structured
        nrow, ncol, nlay, _ = model.nrow_ncol_nlay_nper

        self.itrnsp = itrnsp
        self.ibctcb = ibctcb
        self.mcomp = mcomp
        self.ic_ibound_flg = ic_ibound_flg
        self.itvd = itvd
        self.iadsorb = iadsorb
        self.ict = ict
        self.cinact = cinact
        self.ciclose = ciclose
        self.idisp = idisp
        self.ixdisp = ixdisp
        self.diffnc = diffnc
        self.izod = izod
        self.ifod = ifod
        self.ifmbc = ifmbc
        self.iheat = iheat
        self.imcomp = imcomp
        self.idispcln = idispcln
        self.nseqitr = nseqitr
        self.options = options
        self.icbund = Util3d(
            model,
            (nlay, nrow, ncol),
            np.float32,
            icbund,
            "icbund",
        )
        self.porosity = Util3d(
            model,
            (nlay, nrow, ncol),
            np.float32,
            porosity,
            "porosity",
        )
        self.bulkd = Util3d(
            model,
            (nlay, nrow, ncol),
            np.float32,
            bulkd,
            "bulkd",
        )
        # dispersion settings retrieved from dispd
        self.dispd = dispd
        self.anglex = None
        self.dl = None
        self.dt = None
        self.dlx = None
        self.dly = None
        self.dlz = None
        self.dtxy = None
        self.dlyz = None
        self.dlxz = None
        if self.dispd and idisp != 0:
            if not structured:
                self.anglex = Util2d(
                    model,
                    (dis.njag,),
                    np.float32,
                    dispd["anglex"],
                    "anglex",
                )
            if idisp == 1:
                self.dl = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dl"],
                    "dl",
                )
                self.dt = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dt"],
                    "dt",
                )
            elif idisp == 2:
                self.dlx = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dlx"],
                    "dlx",
                )
                self.dly = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dly"],
                    "dly",
                )
                self.dlz = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dlz"],
                    "dlz",
                )
                self.dtxy = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dtxy"],
                    "dtxy",
                )
                self.dtyz = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dtyz"],
                    "dtyz",
                )
                self.dtxz = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    dispd["dtxz"],
                    "dtxz",
                )

        if self.iadsorb != 0:
            self.adsorb = Util3d(
                model,
                (nlay, nrow, ncol),
                np.float32,
                adsorb,
                "adsorb",
            )
        else:
            self.adsorb = None
        self.sconc = Util3d(
            model,
            (nlay, nrow, ncol),
            np.float32,
            sconc,
            "sconc",
        )
        self.awadsorbd = awadsorbd
        self.awamax = None
        self.alang = None
        self.blang = None
        self.sigma_rt = None
        aw_adsorb = [x.lower() in "a-w_adsorb" for x in self.options]
        if self.awadsorbd and any(aw_adsorb):
            aw_adsorb_i = [i for i, x in enumerate(aw_adsorb) if x]
            iarea_fn = int(self.options[aw_adsorb_i[0] + 1])
            ikawi_fn = int(self.options[aw_adsorb_i[0] + 2])
            if iarea_fn == 1:
                self.awamax = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awadsorbd["awamax"],
                    "awamax",
                )
            elif iarea_fn == 4:
                self.awarea_x2 = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awadsorbd["awarea_x2"],
                    "awarea_x2",
                )
                self.awarea_x1 = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awadsorbd["awarea_x1"],
                    "awarea_x1",
                )
                self.awarea_x0 = Util3d(
                    model,
                    (nlay, nrow, ncol),
                    np.float32,
                    awadsorbd["awarea_x0"],
                    "awarea_x0",
                )
            self.alang = Util3d(
                model,
                (nlay, nrow, ncol),
                np.float32,
                awadsorbd["alang"],
                "alang",
            )
            self.blang = Util3d(
                model,
                (nlay, nrow, ncol),
                np.float32,
                awadsorbd["blang"],
                "blang",
            )
            if ikawi_fn == 3:
                self.sigma_rt = float(awadsorbd["sigma_rt"])

        self.unitnumber_flowtransport = unitnumber_flowtransport
        self.parent.add_package(self)

    def write_file(self, f=None):
        """
        Write the package file.

        Returns
        -------
        None

        """
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper

        # Open file for writing
        if f is None:
            fow = open(self.fn_path, "w")
        else:
            fow = f

        fow.write(f"{self.heading}\n")

        # Item 1: ITRNSP, IBCTCB, MCOMP, IC_IBOUND_FLG, ITVD, IADSORB,
        #         ICT, CINACT, CICLOSE, IDISP, IXDISP, DIFFNC, IZOD, IFOD
        line_1a_vars = (
            self.itrnsp,
            self.ibctcb,
            self.mcomp,
            self.ic_ibound_flg,
            self.itvd,
            self.iadsorb,
            self.ict,
            self.cinact,
            self.ciclose,
            self.idisp,
            self.ixdisp,
            self.diffnc,
            self.izod,
            self.ifod,
            self.ifmbc,
            self.iheat,
            self.imcomp,
            self.idispcln,
            self.nseqitr,
        )
        if not self.parent.free_format_input:
            fstr_format = " <9"
        else:
            fstr_format = ""

        line_1a = [f"{i:{fstr_format}}" for i in line_1a_vars]
        if self.options is not None:
            line_1a += self.options

        fow.write(" ".join(line_1a) + "\n")

        if self.ifmbc != 0:
            line_1b_vars = (
                self.mbegwurf,
                self.mbegwunt,
                self.mbeclnunf,
                self.mbeclnunt,
            )
            line_1b = [f"{i:{fstr_format}}" for i in line_1b_vars]
            fow.write(" ".join(line_1b) + "\n")

        if self.sigma_rt is not None:
            fow.write(f"{self.sigma_rt}" + "\n")  # line 1g (if ikawi_fn == 3)

        if self.ic_ibound_flg == 0:  # line 2
            [
                fow.write(self.icbund[lay].get_file_entry())
                for lay in range(nlay)
            ]

        # line 3
        [fow.write(self.porosity[lay].get_file_entry()) for lay in range(nlay)]

        if self.iadsorb != 0 or self.iheat != 0:  # line 4
            [
                fow.write(self.bulkd[lay].get_file_entry())
                for lay in range(nlay)
            ]

        if self.idisp != 0 and not self.parent.structured:  # line 5
            fow.write(self.anglex.get_file_entry())
            if self.idisp == 1:
                [
                    fow.write(self.dl[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dt[lay].get_file_entry())
                    for lay in range(nlay)
                ]
            elif self.idisp == 2:
                [
                    fow.write(self.dlx[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dly[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dlz[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dtxy[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dtyz[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.dtxz[lay].get_file_entry())
                    for lay in range(nlay)
                ]

        aw_adsorb = [x.lower() in "a-w_adsorb" for x in self.options]
        if any(aw_adsorb):
            aw_adsorb_i = [i for i, x in enumerate(aw_adsorb) if x]
            iarea_fn = int(self.options[aw_adsorb_i[0] + 1])
            if iarea_fn == 1:
                [
                    fow.write(self.awamax[lay].get_file_entry())
                    for lay in range(nlay)
                ]
            elif iarea_fn == 4:
                [
                    fow.write(self.awarea_x2[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.awarea_x1[lay].get_file_entry())
                    for lay in range(nlay)
                ]
                [
                    fow.write(self.awarea_x0[lay].get_file_entry())
                    for lay in range(nlay)
                ]
            [
                fow.write(self.alang[lay].get_file_entry())
                for lay in range(nlay)
            ]
            [
                fow.write(self.blang[lay].get_file_entry())
                for lay in range(nlay)
            ]

        if self.iadsorb != 0:
            [
                fow.write(self.adsorb[lay].get_file_entry())
                for lay in range(nlay)
            ]

        [fow.write(self.sconc[lay].get_file_entry()) for lay in range(nlay)]

        fow.close()

    @classmethod
    def load(cls, f, model, ext_unit_dict=None, check=True):
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading bct package file...")

        fo = get_open_file_object(f, "r")
        eud = ext_unit_dict

        dis = model.get_package("DIS")
        if dis is None:
            dis = model.get_package("DISU")
            nlay = model.nlay
        else:
            nrow, ncol, nlay, _ = dis.nrow_ncol_nlay_nper

        # dataset 0 -- header
        while True:
            line = fo.readline()
            if line[0] != "#":
                break

        if model.verbose:
            print("   loading ITRNSP, IBCTCB, MCOMP...")
        # Item 1a: ITRNSP, IBCTCB, MCOMP, ICBNDFLG, ITVD, IADSORB, ICT, CINACT, CICLOSE, IDISP, IXDISP, DIFFNC, IZOD, IFOD, IFMBC
        text_list = line_parse(line)
        (
            itrnsp,
            ibctcb,
            mcomp,
            ic_ibound_flg,
            itvd,
            iadsorb,
            ict,
            cinact,
            ciclose,
            idisp,
            ixdisp,
            diffnc,
            izod,
            ifod,
            ifmbc,
        ) = (
            int(text_list[0]),
            int(text_list[1]),
            int(text_list[2]),
            int(text_list[3]),
            int(text_list[4]),
            int(text_list[5]),
            int(text_list[6]),
            float(text_list[7]),
            float(text_list[8]),
            int(text_list[9]),
            int(text_list[10]),
            float(text_list[11]),
            int(text_list[12]),
            int(text_list[13]),
            int(text_list[14]),
        )
        iheat, imcomp, idispcln, nseqitr = (None, None, None, None)
        if len(text_list) > 16:
            iheat, imcomp, idispcln, nseqitr = (
                int(x) for x in text_list[15:19]
            )
        options = text_list[19:] if len(text_list) > 19 else None
        aw_adsorb = [False]
        ikawi_fn = None
        if options is not None:
            aw_adsorb = [x.lower() in "a-w_adsorb" for x in options]
            aw_adsorb_i = [i for i, x in enumerate(aw_adsorb) if x]
            if any(aw_adsorb):
                iarea_fn = int(options[aw_adsorb_i[0] + 1])
                ikawi_fn = int(options[aw_adsorb_i[1] + 1])

        # 1b.
        if ifmbc != 0:
            text_list = line_parse(line)
            mbegwurf, mbegwunt, mbeclnunf, mbeclnunt = (
                int(x) for x in text_list
            )
        else:
            mbegwurf, mbegwunt, mbeclnunf, mbeclnunt = (None, None, None, None)
        unitnumber_flowtransport = (mbegwurf, mbegwunt, mbeclnunf, mbeclnunt)

        if ikawi_fn == 3:
            sigma_rt = float(line_parse(line)[0])

        # 2 or 21
        icbund = 0
        if ic_ibound_flg == 0:
            if model.verbose:
                print("   loading icbund...")
            icbund = [
                Util2d.load(
                    fo, model, get_u2d_shape(model, lay), np.int, "icbund", eud
                )
                for lay in range(nlay)
            ]

        # 3 or 22
        if model.verbose:
            print("   loading porosity...")
        porosity = [
            Util2d.load(
                fo,
                model,
                get_u2d_shape(model, lay),
                np.float32,
                "porosity",
                eud,
            )
            for lay in range(nlay)
        ]

        # 4 or 23
        if iadsorb != 0 or (iheat != 0 and not model.structured):
            if model.verbose:
                print("   loading bulkd...")
            bulkd = [
                Util2d.load(
                    fo,
                    model,
                    get_u2d_shape(model, lay),
                    np.float32,
                    "bulkd",
                    eud,
                )
                for lay in range(nlay)
            ]

        # 5
        dispd = {}
        if not model.structured and idisp != 0:
            if model.verbose:
                print("   loading anglex...")
            anglex = Util2d.load(
                fo,
                model,
                (dis.njag,),
                np.float32,
                "anglex",
                eud,
            )
            dispd["anglex"] = anglex
            if idisp == 1:
                # 6 or 24
                if model.verbose:
                    print("   loading dl...")
                dl = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dl",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dl"] = dl

                # 7 or 25
                if model.verbose:
                    print("   loading dt...")
                dt = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dt",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dt"] = dt

            elif idisp == 2:
                # 8 or 26
                if model.verbose:
                    print("   loading dlx...")
                dlx = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dlx",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dlx"] = dlx
                # 9 or 27
                if model.verbose:
                    print("   loading dly...")
                dly = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dly",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dly"] = dly
                # 10 or 28
                if model.verbose:
                    print("   loading dlz...")
                dlz = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dlz",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dlz"] = dlz
                # 11 or 29
                if model.verbose:
                    print("   loading dtxy...")
                dtxy = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dtxy",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dtxy"] = dtxy
                # 12 or 30
                if model.verbose:
                    print("   loading dtyz...")
                dtyz = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dtyz",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dtyz"] = dtyz
                # 13 or 31
                if model.verbose:
                    print("   loading dtyz...")
                dtxz = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "dtxz",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                dispd["dtxz"] = dtxz

        awadsorbd = {}
        if any(aw_adsorb):
            if iarea_fn == 1:
                awamax = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "awamax",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                awadsorbd["awamax"] = awamax
            elif iarea_fn == 4:
                awarea_x2 = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "awarea_x2",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                awadsorbd["awarea_x2"] = awarea_x2
                awarea_x1 = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "awarea_x1",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                awadsorbd["awarea_x1"] = awarea_x1
                awarea_x0 = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "awarea_x0",
                        eud,
                    )
                    for lay in range(nlay)
                ]
                awadsorbd["awarea_x0"] = awarea_x0
            alang = [
                Util2d.load(
                    fo,
                    model,
                    get_u2d_shape(model, lay),
                    np.float32,
                    "alang",
                    eud,
                )
                for lay in range(nlay)
            ]
            awadsorbd["alang"] = alang
            blang = [
                Util2d.load(
                    fo,
                    model,
                    get_u2d_shape(model, lay),
                    np.float32,
                    "blang",
                    eud,
                )
                for lay in range(nlay)
            ]
            awadsorbd["blang"] = blang
            if ikawi_fn == 3:
                awadsorbd["sigma_rt"] = sigma_rt

        for icomp in range(mcomp):
            if icomp == 1:
                print("can only handle 1 species")
                break
            # 14 or 32
            if iadsorb != 0:
                if model.verbose:
                    print(f"   loading adsorb for species {icomp}...")
                adsorb = [
                    Util2d.load(
                        fo,
                        model,
                        get_u2d_shape(model, lay),
                        np.float32,
                        "adsorb",
                        eud,
                    )
                    for lay in range(nlay)
                ]

            if model.verbose:
                print(f"   loading conc for species {icomp}...")
            conc = [
                Util2d.load(
                    fo,
                    model,
                    get_u2d_shape(model, lay),
                    np.float32,
                    "conc",
                    eud,
                )
                for lay in range(nlay)
            ]

        bct = cls(
            model=model,
            itrnsp=itrnsp,
            ibctcb=ibctcb,
            mcomp=mcomp,
            ic_ibound_flg=ic_ibound_flg,
            itvd=itvd,
            iadsorb=iadsorb,
            ict=ict,
            cinact=cinact,
            ciclose=ciclose,
            idisp=idisp,
            ixdisp=ixdisp,
            diffnc=diffnc,
            izod=izod,
            ifod=ifod,
            ifmbc=ifmbc,
            iheat=iheat,
            imcomp=imcomp,
            idispcln=idispcln,
            nseqitr=nseqitr,
            options=options,
            icbund=icbund,
            porosity=porosity,
            bulkd=bulkd,
            dispd=dispd,
            adsorb=adsorb,
            sconc=conc,
            awadsorbd=awadsorbd,
            extension="bct",
            unitnumber=None,
            unitnumber_flowtransport=unitnumber_flowtransport,
        )

        return bct

    @staticmethod
    def _ftype():
        return "BCT"

    @staticmethod
    def _defaultunit():
        return 45
