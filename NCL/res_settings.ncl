
; flag for 2D graph
CUT_BORDER = 14  ; 0/14
LANDSCAPE  = 1  ; 0/1: 0 = portrait, 1 = landscape (better for PPT)

VAR_TYPE   = "PR_MEAN"

COLOR_BAR  = "gui_default"
COLOR_BAR  = "posneg_1"
COLOR_BAR  = "BlWhRe"

LABEL_BAR_ON = True

; var for 1D graph
LINE_THICKNESS = 2.0
LINE_COLOR     = (/"black", "red", "blue", "green4",  "purple", "seagreen4", "brown", "orange", "gold"/)

; 8 region
REG_NUM  = 8
REG_NAME = (/"Cascade", "North Rockies", "Central Great Plain", "Midwest"  \
    ,         "Northeast", "Southeast", "Coast States", "NAM"/)
REG_WDW  = (/(/ 20,  85,  36, 121/)  \  ; (xb, yb, xe, ye) 
    ,        (/ 42,  89,  62, 117/)  \
    ,        (/ 74,  55,  91,  90/)  \
    ,        (/104,  66, 133,  93/)  \
    ,        (/139,  80, 166, 116/)  \
    ,        (/134,  28, 155,  69/)  \
    ,        (/ 94,  39, 128,  63/)  \
    ,        (/ 51,  16,  70,  70/)/)
REG_WDW  = REG_WDW - 1


; output file format
FILE_FORMAT = "pdf"  ; ps, pdf, x11


; shared resource for both 1D and 2D
res_1d = True
res_1d@gsnDraw  = False
res_1d@gsnFrame = False
res_1d@txFont   = "Helvetica"
; res_1d@vpHeightF = 0.4
; res_1d@vpWidthF  = 0.5

res_2d = True
res_2d@gsnDraw  = False
res_2d@gsnFrame = False
res_2d@txFont   = "Helvetica"
; res_2d@vpHeightF = 0.4
; res_2d@vpWidthF  = 0.8


; string and its font height for both 1D and 2D
res_1d@gsnStringFontHeightF       = 0.02
res_1d@gsnLeftStringFontHeightF   = 0.02
res_1d@gsnCenterStringFontHeightF = 0.02
res_1d@gsnRightStringFontHeightF  = 0.02
res_1d@tiMainFontHeightF          = 0.02
res_1d@tiXAxisFontHeightF         = 0.02
res_1d@tiYAxisFontHeightF         = 0.02

res_1d@gsnLeftString   = ""
res_1d@gsnCenterString = ""
res_1d@gsnRightString  = ""
res_1d@tiMainString    = ""
res_1d@tiXAxisString   = ""
res_1d@tiYAxisString   = ""
res_1d@tiXAxisSide     = "Bottom" 
res_1d@tiYAxisSide     = "Left"

res_2d@gsnStringFontHeightF       = 0.02
res_2d@gsnLeftStringFontHeightF   = 0.02
res_2d@gsnCenterStringFontHeightF = 0.02
res_2d@gsnRightStringFontHeightF  = 0.02
res_2d@tiMainFontHeightF          = 0.02
res_2d@tiXAxisFontHeightF         = 0.02
res_2d@tiYAxisFontHeightF         = 0.02

res_2d@gsnLeftString   = ""
res_2d@gsnCenterString = ""
res_2d@gsnRightString  = ""
res_2d@tiMainString    = ""
res_2d@tiXAxisString   = ""
res_2d@tiYAxisString   = ""
res_2d@tiXAxisSide     = "Bottom"
res_2d@tiYAxisSide     = "Left"


; tickmark for both 1D and 2D
res_1d@pmTickMarkDisplayMode   = "Always"  ; Never, NoCreate, Conditional
res_1d@tmXTOn                  = False
res_1d@tmYROn                  = False
res_1d@tmXBMinorOn             = False
res_1d@tmYLMinorOn             = False

res_1d@tmXBMajorLengthF        = 0.005
res_1d@tmYLMajorLengthF        = 0.005
res_1d@tmXBMajorOutwardLengthF = -0.002  ; draw tickmarks inward
res_1d@tmYLMajorOutwardLengthF = -0.002
res_1d@tmXBLabelDeltaF         = -1.2
res_1d@tmYLLabelDeltaF         = -1.2
res_1d@tmXBLabelFontHeightF    = 0.01
res_1d@tmYLLabelFontHeightF    = 0.01

res_1d@tmXBMode     = "Explicit"
res_1d@tmXBLabelsOn = True
res_1d@tmXBValues   = ispan(0, 11, 1)
res_1d@tmXBLabels   = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"/)

; res_1d@tmYLMode          = "Manual"
; res_1d@tmYLLabelsOn      = True
; res_1d@trYMinF           = 0.0
; res_1d@trYMaxF           = 18.0
; res_1d@tmYLTickStartF    = 0.0
; res_1d@tmYLTickEndF      = 18.0
; res_1d@tmYLTickSpacingF  = 6.0 
; res_1d@tmYLAutoPrecision = False
; res_1d@tmYLPrecision     = 2

res_2d@pmTickMarkDisplayMode   = "Never"  ; "Always", NoCreate, Conditional
res_2d@tmXTOn                  = False
res_2d@tmYROn                  = False
res_2d@tmXBMinorOn             = False
res_2d@tmYLMinorOn             = False

res_2d@tmXBMajorLengthF        = 0.005
res_2d@tmYLMajorLengthF        = 0.005
res_2d@tmXBMajorOutwardLengthF = -0.002  ; draw tickmarks inward
res_2d@tmYLMajorOutwardLengthF = -0.002
res_2d@tmXBLabelDeltaF         = -1.2
res_2d@tmYLLabelDeltaF         = -1.2
res_2d@tmXBLabelFontHeightF    = 0.01
res_2d@tmYLLabelFontHeightF    = 0.01

; res_2d@tmXBMode     = "Explicit"
; res_2d@tmXBLabelsOn = True
; res_2d@tmXBValues   = ispan(0, 11, 1)
; res_2d@tmXBLabels   = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"/)

; res_2d@tmYLMode          = "Manual"
; res_2d@tmYLLabelsOn      = True
; res_2d@trYMinF           = 0.0
; res_2d@trYMaxF           = 18.0
; res_2d@tmYLTickStartF    = 0.0
; res_2d@tmYLTickEndF      = 18.0
; res_2d@tmYLTickSpacingF  = 6.0 
; res_2d@tmYLAutoPrecision = False
; res_2d@tmYLPrecision     = 2

; res_2d@tmYLMode       = "Explicit"  ; for ETS plot
; res_2d@tmYLValues     = (/5  , 11 , 17 , 23  , 29  /)
; res_2d@tmYLLabels     = (/"3", "6", "9", "12", "15"/)


; legend and line/bar colors and patterns for 1D
res_1d@pmLegendDisplayMode    = "Never"  ; "Always" 
res_1d@pmLegendSide           = "Top"  ; "Bottom"
res_1d@pmLegendParallelPosF   = 0.6    ; bigger = right
res_1d@pmLegendOrthogonalPosF = -0.6   ; more negative = down
res_1d@pmLegendWidthF         = 0.1 
res_1d@pmLegendHeightF        = 0.15

res_1d@lgLabelFontHeightF     = 0.015
res_1d@lgPerimOn              = False
; res_1d@lgItemOrder            = (/4, 3, 2, 1, 0/)

; res_1d@xyExplicitLegendLabels = (/"CWRF/ERI", "CWRF", "WRF", "WRFG", "CMM5"/)
; res_1d@xyLineThicknesses      = (/2.0, 2.0, 2.0, 2.0, 2.0/)
; res_1d@xyDashPatterns         = (/1, 0, 0, 1, 0/)
; res_1d@xyLineColors           = (/"red", "red", "(/0.498, 1.0, 0.831/)", "cyan2", "blue"/)

res_1d@xyExplicitLegendLabels = (/"AVG", "OPT"/)
res_1d@xyLineThicknesses      = (/2.0, 2.0/)
res_1d@xyDashPatterns         = (/0, 0/)
res_1d@xyLineColors           = (/"blue", "red"/)

; color maps for 2D
; gsn_define_colormap(wks_pdf, "gui_default")  

; gsn_merge_colormaps(wks_pdf, "gui_default", "posneg_1")

; colors = (/"cyan2", "(/0.498, 1.00, 0.831/)"/)
; gsn_merge_colormaps(wks_pdf, "gui_default", colors)  ; necessary to use cyan2

; gsn_define_colormap(wks_pdf, "gui_default")
; colors      = gsn_retrieve_colormap(wks_pdf)
; cmap        = new ((/12, 3/), "float")
; cmap(0 , :) = colors(0 , :)
; cmap(1 , :) = colors(1 , :)
; cmap(2 , :) = colors(3 , :)
; cmap(3 , :) = colors(5 , :)
; cmap(4 , :) = colors(9 , :)
; cmap(5 , :) = colors(17, :)
; cmap(6 , :) = colors(20, :)
; cmap(7 , :) = colors(21, :)
; cmap(8 , :) = colors(22, :)
; cmap(9 , :) = colors(23, :)
; cmap(10, :) = (/125.0, 38.0, 205.0/) / 255.0  ; purple3
; cmap(11, :) = (/85.0 , 26.0, 139.0/) / 255.0  ; purple4
; gsn_define_colormap(wks_pdf, cmap)


; maps for 2D
res_2d@mpLimitMode           = "Corners"
if(CUT_BORDER .eq. 0) then
    res_2d@mpLeftCornerLatF  =   15.22807
    res_2d@mpLeftCornerLonF  = -120.7703
    res_2d@mpRightCornerLatF =   49.03231
    res_2d@mpRightCornerLonF =  -52.39819
else if(CUT_BORDER .eq. 14) then
    res_2d@mpLeftCornerLatF  =   19.768133
    res_2d@mpLeftCornerLonF  = -118.259
    res_2d@mpRightCornerLatF =   47.4818
    res_2d@mpRightCornerLonF =  -60.0576
else
      print("ERROR: Unknow CUT_BORDER VALUE, PLEASE CHECK.")
      exit
end if
end if

res_2d@mpProjection          = "LambertConformal"
res_2d@mpLambertParallel1F   = 60.
res_2d@mpLambertParallel2F   = 30.
res_2d@mpLambertMeridianF    = -95.5

res_2d@mpFillOn              = False
res_2d@mpOutlineDrawOrder    = "PostDraw"
res_2d@mpOutlineBoundarySets = "GeophysicalAndUSStates"

res_2d@gsnAddCyclic             = False
res_2d@tfDoNDCOverlay           = True  ; if False, then "check_for_y_lat_coord: Warning: Data either ..."


; contour and label bar for 2D
res_2d@cnLinesOn      = False
res_2d@cnLineLabelsOn = False
res_2d@cnFillOn       = True
res_2d@cnInfoLabelOn  = False

res_2d@gsnSpreadColors     = True
; res_2d@gsnSpreadColorStart = 2
; res_2d@gsnSpreadColorEnd   = 23
; res_2d@gsnSpreadColorStart = 27
; res_2d@gsnSpreadColorEnd   = 44

res_2d@cnLevelSelectionMode = "ManualLevels"
if (VAR_TYPE .eq. "PR_MEAN") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 15.0
    res_2d@cnLevelSpacingF = 1.0
else if (VAR_TYPE .eq. "PR_MEAN_DIFF") then
    res_2d@cnMinLevelValF  = -5.0
    res_2d@cnMaxLevelValF  = 5.0
    res_2d@cnLevelSpacingF = 0.5 
else if (VAR_TYPE .eq. "T2M_MEAN") then
    res_2d@cnMinLevelValF  = -10.0
    res_2d@cnMaxLevelValF  = 35.0
    res_2d@cnLevelSpacingF = 1.0
else if (VAR_TYPE .eq. "T2M_MEAN_DIFF") then
    res_2d@cnMinLevelValF  = -5.0
    res_2d@cnMaxLevelValF  = 5.0
    res_2d@cnLevelSpacingF = 0.5
else if (VAR_TYPE .eq. "SSN_CC") then
    res_2d@cnMinLevelValF  = -0.2
    res_2d@cnMaxLevelValF  = 1.0
    res_2d@cnLevelSpacingF = 0.1
else if (VAR_TYPE .eq. "SSN_CC_DIFF") then
    res_2d@cnMinLevelValF  = -0.3
    res_2d@cnMaxLevelValF  = 0.3
    res_2d@cnLevelSpacingF = 0.05
else if (VAR_TPYE .eq. "SSN_RMS") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 15.0
    res_2d@cnLevelSpacingF = 1.0 
else if (VAR_TPYE .eq. "SSN_RMS_DIFF") then
    res_2d@cnMinLevelValF  = -5.0
    res_2d@cnMaxLevelValF  = 5.0
    res_2d@cnLevelSpacingF = 1.0
else if (VAR_TYPE .eq. "ETS") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 1.0
    res_2d@cnLevelSpacingF = 0.05
else if (VAR_TYPE .eq. "ETS_DIFF") then
    res_2d@cnMinLevelValF  = -0.4
    res_2d@cnMaxLevelValF  = 0.4
    res_2d@cnLevelSpacingF = 0.05
else if (VAR_TYPE .eq. "RPS") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 0.2
    res_2d@cnLevelSpacingF = 0.01
else if (VAR_TYPE .eq. "RPS_DIFF") then
    res_2d@cnMinLevelValF  = -0.02
    res_2d@cnMaxLevelValF  = 0.02
    res_2d@cnLevelSpacingF = 0.002
else if (VAR_TYPE .eq. "RPSS") then
    res_2d@cnMinLevelValF  = -0.1
    res_2d@cnMaxLevelValF  = 1.0 
    res_2d@cnLevelSpacingF = 0.05 
else if (VAR_TYPE .eq. "RPSS_DIFF") then
    res_2d@cnMinLevelValF  = -0.3
    res_2d@cnMaxLevelValF  = 0.3 
    res_2d@cnLevelSpacingF = 0.05
else if (VAR_TYPE .eq. "SSN_NUM_OF_RAINY_DAYS") then
    res_2d@cnMinLevelValF  = 8.0 
    res_2d@cnMaxLevelValF  = 40.0
    res_2d@cnLevelSpacingF = 2.0 
else if (VAR_TYPE .eq. "SSN_MAX_DRY_SPELL_LENGTH") then
    res_2d@cnMinLevelValF  = 10.0
    res_2d@cnMaxLevelValF  = 70.0
    res_2d@cnLevelSpacingF = 5.0
else if (VAR_TYPE .eq. "SSN_DAILY_RAINFALL_95PERCENTILE") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 30.0
    res_2d@cnLevelSpacingF = 2.0
else if (VAR_TYPE .eq. "SSN_PBLH") then
    res_2d@cnMinLevelValF  = 0.0
    res_2d@cnMaxLevelValF  = 2000.0
    res_2d@cnLevelSpacingF = 100.0
else if (VAR_TYPE .eq. "SSN_PBLH_DIFF") then
    res_2d@cnMinLevelValF  = -500.0
    res_2d@cnMaxLevelValF  = 500.0
    res_2d@cnLevelSpacingF = 50.0
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if

if(LABEL_BAR_ON) then
   res_2d@lbLabelBarOn = True 
else
   res_2d@lbLabelBarOn = False
end if
res_2d@lbOrientation            = "Horizontal"  ; "Vertical"
res_2d@pmLabelBarOrthogonalPosF = -0.01  ; move label bar closer
res_2d@lbLabelAutoStride        = True 
res_2d@lbLabelStride            = 2
res_2d@lbAutoManage             = False
res_2d@lbLabelFontHeightF       = 0.015


; panel resource for function gsn_panel
res_panel               = True
res_panel@gsnPanelDebug = True  ; To get information from panel

; res_panel@gsnPanelXF = (/0.194416, 0.496775, 0.194416, 0.496775, 0.194416, 0.496775, 0.194416, 0.496775/)
; res_panel@gsnPanelYF = (/0.966156, 0.966156, 0.767619, 0.767619, 0.569082, 0.569082, 0.370545, 0.370545/)
; res_panel@gsnPanelXF = (/0.123753, 0.365764, 0.607775, 0.123753, 0.365764, 0.607775, 0.123753, 0.365764, 0.607775  \
;     , 0.123753, 0.365764, 0.607775, 0.123753, 0.365764, 0.607775/)
; res_panel@gsnPanelYF = (/0.973013, 0.973013, 0.973013, 0.814102, 0.814102, 0.814102, 0.655191, 0.655191, 0.655191  \
;     , 0.496280, 0.496280, 0.496280, 0.337369, 0.337369, 0.337369/)

res_panel@gsnPanelLabelBar         = False
res_panel@lbOrientation            = "Horizontal"
res_panel@pmLabelBarOrthogonalPosF = -0.01  ; 0.0, the bigger, the closer of the label bar
res_panel@lbLabelAutoStride        = True   ; optimal labels
res_panel@lbLabelStride            = 2
res_panel@lbAutoManage             = False
res_panel@lbLabelFontHeightF       = 0.015

; res_panel@cnLevelSelectionMode   = "ExplicitLevels"
; res_panel@cnLevels               = fspan(-4.0, 4.0, 16)

; res_panel@cnLevelSelectionMode   = "ManualLevels"  ; wouldn't work!
; res_panel@cnLevelSpacingF        = 0.5
; res_panel@cnMinLevelValF         = 0.0
; res_panel@cnMaxLevelValF         = 12.0


; text resource for function gsn_add_text to add text to plot, and an example
res_text               = True
res_text@txFontHeightF = 0.015
res_text@txAngleF      = 0.0
res_text@txDirection   = "Across"
res_text@txFontColor   = 1

; getvalues plot_pdf(0)
;     "tmXTValues": tmXTValues
;     "tmYLValues": tmYLValues
; end getvalues
; nXm = dimsizes(tmXTValues)
; nYm = dimsizes(tmYLValues)

; TEXT = (/":F21:PR:N: :C: :F34:/:N:", " :F34:-:N: :C:" + ":F21:T2M:N:"/)
; d    = gsn_add_text(wks_pdf, plot_pdf(0), TEXT, (/0.65*tmXTValues(nXm-1), 0.45*tmXTValues(nXm-1)/)  \
;     , (/0.65*tmYLValues(nYm-1), 0.85*tmYLValues(nYm-1)/), res_text)


; polygon resource for function gsn_add_polygon to add polygon to plot, and an example
res_polygon                   = True
res_polygon@gsFillColor       = -1
res_polygon@gsEdgesOn         = True
res_polygon@gsEdgeColor       = "black"
res_polygon@gsEdgeThicknessF  = 2.0
res_polygon@gsEdgeDashPattern = 0

; x_pg = (/lon(yb, xb), lon(yb, xe), lon(ye, xe), lon(ye, xb), lon(yb, xb)/)
; y_pg = (/lat(yb, xb), lat(yb, xe), lat(ye, xe), lat(ye, xb), lat(yb, xb)/)
; pg(0, i_reg) = gsn_add_polygon(wks_pdf, plot_pdf(0), x_pg, y_pg, res_polygon)


