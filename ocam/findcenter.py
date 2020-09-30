###########################################################################  
#   Copyright (C) 2006 DAVIDE SCARAMUZZA
#   
#   Author: Davide Scaramuzza - email: davsca@tiscali.it
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#   USA
############################################################################

from numpy import meshgrid

from .calib_data import CalibData
from .calibrate import calibrate
from .libsmop import *
from .reprojectpoints import reprojectpoints
from .reprojectpoints_fun import reprojectpoints_fun


@function
def findcenter(calib_data: CalibData):
    if isempty(calib_data.ima_proc) or isempty(calib_data.Xp_abs):
        print('\nNo corner data available. Extract grid corners before calibrating.\n')
        return

    print('\nComputing center coordinates.\n')
    if isempty(calib_data.taylor_order):
        calib_data.taylor_order = calib_data.taylor_order_default

    pxc = calib_data.ocam_model.xc
    pyc = calib_data.ocam_model.yc
    width = calib_data.width
    height = calib_data.height
    regwidth = (width / 2)
    regheight = (height / 2)
    yceil = 5
    xceil = 5
    xregstart = pxc - (regheight / 2)
    xregstop = pxc + (regheight / 2)
    yregstart = pyc - (regwidth / 2)
    yregstop = pyc + (regwidth / 2)
    print('Iteration ', end='')
    for glc in range(1, 9 + 1):
        yreg, xreg = meshgrid(arange(yregstart, yregstop + 1 / yceil, (yregstop - yregstart) / yceil),
                              arange(xregstart, xregstop + 1 / xceil, (xregstop - xregstart) / xceil), nargout=2)
        ic_proc = concat([arange(1, size(xreg, 1))])
        jc_proc = concat([arange(1, size(xreg, 2))])
        MSEA = dot(inf, ones(size(xreg)))
        for ic in ic_proc.flat:
            for jc in jc_proc.flat:
                calib_data.ocam_model.xc = xreg(ic, jc)
                calib_data.ocam_model.yc = yreg(ic, jc)
                #            hold on; plot(yc,xc,'r.');
                calib_data.RRfin, calib_data.ocam_model.ss = calibrate(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs,
                                                                       calib_data.Yp_abs, calib_data.ocam_model.xc,
                                                                       calib_data.ocam_model.yc,
                                                                       calib_data.taylor_order, calib_data.ima_proc,
                                                                       nargout=2)
                if calib_data.RRfin == 0:
                    MSEA[ic, jc] = inf
                    continue
                MSE = reprojectpoints_fun(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs, calib_data.Yp_abs,
                                          calib_data.ocam_model.xc, calib_data.ocam_model.yc, calib_data.RRfin,
                                          calib_data.ocam_model.ss, calib_data.ima_proc, calib_data.width,
                                          calib_data.height)
                # obrand_start 
                # speedup removed to compensate for calibration errors
                #            if ic>1 & jc>1
                #                if MSE>MSEA(ic-1,jc)
                #                    jc_proc(find(jc_proc==jc))=inf;
                #                    jc_proc=sort(jc_proc);
                #                    jc_proc=jc_proc[1:end-1];
                #                    continue;
                #                elseif MSE>MSEA(ic,jc-1)
                #                    break;
                #                elseif isnan(MSE)
                #                    break;
                #                end
                #            end
                #            MSEA(ic,jc)=MSE;
                # obrand_replacement
                if logical_not(isnan(MSE)):
                    MSEA[ic, jc] = MSE
                # obrand_end
        #    drawnow;
        indMSE = find(min(ravel(MSEA)) == MSEA)
        calib_data.ocam_model.xc = xreg(indMSE[1])
        calib_data.ocam_model.yc = yreg(indMSE[1])
        dx_reg = abs((xregstop - xregstart) / xceil)
        dy_reg = abs((yregstop - yregstart) / yceil)
        xregstart = calib_data.ocam_model.xc - dx_reg
        xregstop = calib_data.ocam_model.xc + dx_reg
        yregstart = calib_data.ocam_model.yc - dy_reg
        yregstop = calib_data.ocam_model.yc + dy_reg
        print('%d...', glc, end='')

    print('')
    calib_data.RRfin, calib_data.ocam_model.ss = calibrate(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs,
                                                           calib_data.Yp_abs, calib_data.ocam_model.xc,
                                                           calib_data.ocam_model.yc, calib_data.taylor_order,
                                                           calib_data.ima_proc, nargout=2)
    reprojectpoints(calib_data)
    xc = calib_data.ocam_model.xc
    yc = calib_data.ocam_model.yc
    xc
    yc
    calib_data.calibrated = 1
