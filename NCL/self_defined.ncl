; Description: A bunch of self-defined NCL functions and procedures
; Author     : Ligang Chen, lchen2@umd.edu, ESSIC/UMD


; the result index is origined at 0, which is different from fortran language
function get_record_index(origin_year, origin_month, origin_day, sample_year, sample_month, sample_day)
begin
    dofy_origin = day_of_year(origin_year, origin_month, origin_day)
    dofy_sample = day_of_year(sample_year, sample_month, sample_day)

    res_index = 0

    if (sample_year .eq. origin_year) then
        res_index = dofy_sample - dofy_origin + 1
    else
        if (isleapyear(origin_year)) then
            res_index = res_index + 366 - dofy_origin + 1
        else
            res_index = res_index + 365 - dofy_origin + 1
        end if

        res_index = res_index + dofy_sample

        do i_yr = origin_year + 1, sample_year - 1
            if (isleapyear(i_yr)) then
                res_index = res_index + 366
            else
                res_index = res_index + 365
            end if
        end do
    end if

    res_index = res_index - 1  ; make the index origined at 0

    return res_index
end  ; function get_record_index


; the result index is origined at 0, which is different from fortran language
function get_record_index_no_leap_year(origin_year, origin_month, origin_day, sample_year, sample_month, sample_day)
begin
    dofy_origin = day_of_year(origin_year, origin_month, origin_day)
    dofy_sample = day_of_year(sample_year, sample_month, sample_day)

    if (isleapyear(origin_year) .and. origin_month .ge. 3) then
        dofy_origin = dofy_origin - 1
    end if
    if (isleapyear(sample_year) .and. sample_month .ge. 3) then
        dofy_sample = dofy_sample - 1
    end if

    res_index = 0

    if (sample_year .eq. origin_year) then
        res_index = dofy_sample - dofy_origin + 1
    else
        if (isleapyear(origin_year)) then
            res_index = res_index + 365 - dofy_origin + 1
        else
            res_index = res_index + 365 - dofy_origin + 1
        end if

        res_index = res_index + dofy_sample

        do i_yr = origin_year + 1, sample_year - 1
            if (isleapyear(i_yr)) then
                res_index = res_index + 365
            else
                res_index = res_index + 365
            end if
        end do
    end if

    res_index = res_index - 1  ; make the index origined at 0

    return res_index
end  ; function get_record_index_no_leap_year


; the result index is origined at 0
function get_record_index_monthly(origin_year, origin_month, sample_year, sample_month)
begin
    res_index = 0

    if (sample_year .eq. origin_year) then
        res_index = sample_month - origin_month + 1
    else
        res_index = res_index + 12 - origin_month + 1
        res_index = res_index + sample_month

        do i_yr = origin_year + 1, sample_year - 1
            res_index = res_index + 12
        end do
    end if

    res_index = res_index - 1  ; make the index origined at 0

    return res_index
end  ; function get_record_index_monthly


; the result index is origined at 0
function get_record_index_seasonal(origin_year, origin_season, sample_year, sample_season)
begin
    res_index = 0

    if (sample_year .eq. origin_year) then
        res_index = sample_season - origin_season + 1
    else
        res_index = res_index + 4 - origin_season + 1
        res_index = res_index + sample_season

        do i_yr = origin_year + 1, sample_year - 1
            res_index = res_index + 4
        end do
    end if

    res_index = res_index - 1  ; make the index origined at 0

    return res_index
end  ; function get_record_index_seasonal


; the result index is origined at 0, which is different from fortran language
function get_record_index_6hourly(origin_year, origin_month, origin_day, origin_cycle  \
    , sample_year, sample_month, sample_day, sample_cycle)
begin
    dofy_origin = day_of_year(origin_year, origin_month, origin_day)
    dofy_sample = day_of_year(sample_year, sample_month, sample_day)

    ; res_index = 0
    n_day_interval = 0

    if (sample_year .eq. origin_year) then
        n_day_interval = dofy_sample - dofy_origin + 1
    else
        if (isleapyear(origin_year)) then
            n_day_interval = n_day_interval + 366 - dofy_origin + 1
        else
            n_day_interval = n_day_interval + 365 - dofy_origin + 1
        end if

        n_day_interval = n_day_interval + dofy_sample

        do i_yr = origin_year + 1, sample_year - 1
            if (isleapyear(i_yr)) then
                n_day_interval = n_day_interval + 366
            else
                n_day_interval = n_day_interval + 365
            end if
        end do
    end if

    ; the above n_day_interval includes both the origin_day and the sample_day, exclude them
    n_day_interval = n_day_interval - 2

    if (n_day_interval .lt. 0) then  ; n_day_interval would be equal to -1
        res_index = sample_cycle - origin_cycle
    else
        res_index = 4*n_day_interval + sample_cycle + (5 - origin_cycle) - 1
    end if

    return res_index
end  ; function get_record_index_6hourly


procedure compute_freq_cc_x(freq_cc_x:float, n_seg_cc[1]:integer, x_mid[1]:logical)
; normally n_seg_cc = 40 and bin_len = 0.05
; freq_cc_x(n_seg_cc)
begin
    do i_seg = 1, n_seg_cc
        if (x_mid) then
            freq_cc_x(i_seg - 1) = -1.0 + (2.0*i_seg - 1.0) / n_seg_cc
        else
            freq_cc_x(i_seg - 1) = -1.0 + (2.0*i_seg      ) / n_seg_cc
        end if
    end do
end  ; procedure compute_freq_cc_x


procedure compute_freq_rms_x(freq_rms_x:float, num_step_rms[1]:integer, len_step_rms[1]:float, x_mid[1]:logical)
; normally daily PR     : num_step_rms = 30, len_step_rms = 1.0
; seasonal or monthly PR: num_step_rms = 40, len_step_rms = 0.1
begin
    do i_seg = 1, num_step_rms
        if (x_mid) then
            freq_rms_x(i_seg - 1) = 0.5*len_step_rms + (i_seg-1)*len_step_rms
        else
            freq_rms_x(i_seg - 1) = i_seg*len_step_rms
        end if
    end do
end  ; procedure compute_freq_rms_x


; compute the frequency distribution of one record of CC and RMS
procedure compute_freq_cc_rms(cc[138][195]:float, rms[138][195]:float, freq_cc:float, freq_rms:float  \
    , lcc[138][195]:float, us_landmask[138][195]:float                                      \
    , buffer_zone[1]:logical, land[1]:logical, ocean[1]:logical, us_land_only[1]:logical    \
    , n_seg_cc[1]:integer, len_step_rms[1]:float, num_step_rms[1]:integer  \
    , str_log[500]:string, n_line_log[1]:integer)
; freq_cc(n_seg_cc), freq_rms(num_step_rms)
; normally n_seg_cc = 40, num_step_rms = 30 and len_step_rms = 1.0 for daily PR, 
;     num_step_rms = 40 and len_step_rms = 0.1 for seasonal or monthly mean PR.
begin
    ; initialize
    freq_cc  = 0.0
    freq_rms = 0.0

    n_line_log = 0
    str_log    = ""

    nvp_cc   = 0
    nivp_cc  = 0
    nvp_rms  = 0
    nivp_rms = 0

    do i_y = 0, 138 - 1
        do i_x = 0, 195 - 1
            if (.not. buffer_zone .and. (i_y .lt. 14 .or. i_y .gt. 138 - 1 - 14 .or. i_x .lt. 14  \
                .or. i_x .gt. 195 - 1 -14)) then
                continue
            end if

            if (us_land_only .and. us_landmask(i_y, i_x) .lt. 0.5) then
                continue
            end if

            if (.not. land .and. lcc(i_y, i_x) .ne. 16.0) then
                continue
            end if

            if (.not. ocean .and. lcc(i_y, i_x) .eq. 16.0) then
                continue
            end if

            ; cc
            if (.not. ismissing(cc(i_y, i_x)) .and. -1.0 .le. cc(i_y, i_x) .and. cc(i_y, i_x) .le. 1.0) then
                nvp_cc = nvp_cc + 1

                do i_seg = 1, n_seg_cc
                    if (cc(i_y, i_x) .gt. 2.0*(i_seg-1)/n_seg_cc - 1.0  \
                        .and. cc(i_y, i_x) .le. 2.0*i_seg/n_seg_cc - 1.0) then
                        break
                    end if
                end do

                freq_cc(i_seg - 1) = freq_cc(i_seg - 1) + 1
            else
                nivp_cc = nivp_cc + 1
                str_log(n_line_log) = "Invalid cc: cc(" + sprinti("%0.2i", i_y) + ", " + sprinti("%0.2i", i_x)  \
                    + ") = " + sprintf("%7.2f", cc(i_y, i_x))
                n_line_log = n_line_log + 1
                ; asciiwrite(fn_freq_log, str_log)
            end if

            ; rms
            if (0.0 .le. rms(i_y, i_x) .and. rms(i_y, i_x) .le. 1.0e5) then
                nvp_rms = nvp_rms + 1

                do i_seg = 1, num_step_rms
                    if (1.0*(i_seg-1)*len_step_rms .le. rms(i_y, i_x) .and. rms(i_y, i_x) .lt. 1.0*i_seg*len_step_rms) then
                        break
                    end if
                end do
                if (i_seg .le. num_step_rms) then
                    freq_rms(i_seg - 1) = freq_rms(i_seg - 1) + 1
                else
                    ; freq_rms(num_step_rms - 1) = freq_rms(num_step_rms - 1) + 1  ; would cause sudden leap
                end if
            else
                nivp_rms = nivp_rms + 1
                str_log(n_line_log) = "Invalid rms: rms(" + sprinti("%0.2i", i_y) + ", " + sprinti("%0.2i", i_x)  \
                    + ") = " + sprintf("%7.2f", rms(i_y, i_x))
                n_line_log = n_line_log + 1
            end if
        end do  ; end of "do i_x = 0, 195 - 1"
    end do      ; end of "do i_y = 0, 138 - 1"

    freq_cc  = freq_cc  / (1.0*nvp_cc )
    freq_rms = freq_rms / (1.0*nvp_rms)

    str_log(n_line_log) = "nvp_cc = " + sprinti("%0.6i", nvp_cc) + ", nivp_cc = " + sprinti("%0.6i", nivp_cc)  \
        + ", nvp_rms = " + sprinti("%0.6i", nvp_rms) + ", nivp_rms = " + sprinti("%0.6i", nivp_rms)
    n_line_log = n_line_log + 1

    return
end  ; procedure compute_freq_cc_rms


; compute 1D CC frequency distribution
procedure compute_freq_cc_1d(cc:float, freq_cc:float  \
    , n_cc[1]:integer, n_seg_cc[1]:integer, str_log[500]:string, n_line_log[1]:integer)
; cc(n_cc), freq_cc(n_seg_cc)
begin
    ; initialize
    freq_cc  = 0.0

    n_line_log = 0
    str_log    = ""

    nvp_cc   = 0
    nivp_cc  = 0

    do i_cc = 0, n_cc - 1
        if (.not. ismissing(cc(i_cc)) .and. -1.0 .le. cc(i_cc) .and. cc(i_cc) .le. 1.0) then
            nvp_cc = nvp_cc + 1

            do i_seg = 1, n_seg_cc
                if (cc(i_cc) .gt. 2.0*(i_seg-1)/n_seg_cc - 1.0 .and. cc(i_cc) .le. 2.0*i_seg/n_seg_cc - 1.0) then
                    break
                end if
            end do

            freq_cc(i_seg - 1) = freq_cc(i_seg - 1) + 1
        else
            nivp_cc = nivp_cc + 1
            str_log(n_line_log) = "Invalid cc: cc(" + sprinti("%0.3i", i_cc) + ") = " + sprintf("%7.2f", cc(i_cc))
            n_line_log = n_line_log + 1
        end if
    end do

    freq_cc  = freq_cc  / (1.0*nvp_cc )

    str_log(n_line_log) = "nvp_cc = " + sprinti("%0.6i", nvp_cc) + ", nivp_cc = " + sprinti("%0.6i", nivp_cc)
    n_line_log          = n_line_log + 1

    return
end  ; procedure compute_freq_cc_1d


; compute 1D RMSE frequency distribution
procedure compute_freq_rms_1d(rms:float, freq_rms:float  \
    , n_rms[1]:integer, len_step_rms[1]:float, num_step_rms[1]:integer  \
    , str_log[500]:string, n_line_log[1]:integer)
; rms(n_rms), freq_rms(num_step_rms)
begin
    ; initialize
    freq_rms = 0.0

    n_line_log = 0
    str_log    = ""

    nvp_rms  = 0
    nivp_rms = 0

    do i_rms = 0, n_rms - 1
        if (0.0 .le. rms(i_rms) .and. rms(i_rms) .le. 1.0e5) then
            nvp_rms = nvp_rms + 1

            do i_seg = 1, num_step_rms
                if (1.0*(i_seg-1)*len_step_rms .le. rms(i_rms) .and. rms(i_rms) .lt. 1.0*i_seg*len_step_rms) then
                    break
                end if
            end do
            if (i_seg .le. num_step_rms) then
                freq_rms(i_seg - 1) = freq_rms(i_seg - 1) + 1
            else
                ; freq_rms(num_step_rms - 1) = freq_rms(num_step_rms - 1) + 1  ; would cause sudden leap
            end if
        else
            nivp_rms = nivp_rms + 1
            str_log(n_line_log) = "Invalid rms: rms(" + sprinti("%0.3i", i_rms) + ") = " + sprintf("%7.2f", rms(i_rms))
            n_line_log = n_line_log + 1
        end if
    end do      ; end of "do i_rms = 0, n_rms - 1"

    freq_rms = freq_rms / (1.0*nvp_rms)

    str_log(n_line_log) = "nvp_rms = " + sprinti("%0.6i", nvp_rms) + ", nivp_rms = " + sprinti("%0.6i", nivp_rms)
    n_line_log = n_line_log + 1

    return
end  ; procedure compute_freq_rms_1d


; extend a string to fixed length
function str_fix_len(str, len)
begin
    n_extra_len = len - strlen(str)
    if (n_extra_len .ge. 1) then
        blank_str    = new (n_extra_len, string)
        blank_str(:) = " "
        return (str + str_concat(blank_str))
    else
        return (str)
    end if 
end  ; function str_fix_len


procedure set_8reg_mask(reg8_mask)
; reg8_mask(138, 195)
begin
    N_REG    = 8
    REG_NAME = (/"Cascade", "North Rockies", "Central Great Plain", "Midwest"  \
        , "Northeast", "Southeast", "Coast States", "NAM"/)
    REG_WIN  = (/(/ 20,  85,  36, 121/)  \  ; (xb, yb, xe, ye)
        ,        (/ 42,  89,  62, 117/)  \
        ,        (/ 74,  55,  91,  90/)  \
        ,        (/104,  66, 133,  93/)  \
        ,        (/139,  80, 166, 116/)  \
        ,        (/134,  28, 155,  69/)  \
        ,        (/ 94,  39, 128,  63/)  \
        ,        (/ 51,  16,  70,  70/)/)
    REG_WIN = REG_WIN - 1

    reg8_mask = 0.0

    do i_reg = 0, N_REG - 1
        reg8_mask(REG_WIN(i_reg, 1):REG_WIN(i_reg, 3), REG_WIN(i_reg, 0):REG_WIN(i_reg, 2)) = 1.0
    end do
end  ;  procedure set_8reg_mask


function create_legend_from_plot(wks:graphic, plot:graphic, opt:logical)
local lgres, classname

; This function takes an existing XY plot and constructs a legend from it. 
; This is so you can further customize the legend as desired.
;
; If "opt" is set to True, then it can have two attributes:
; LegendLabels      - an array of strings to use for the legend labels
; NumberLegendItems - the number of items you want to display.
;     The default will be however many legend labels there are.

; Check it out at http://www.ncl.ucar.edu/Applications/legend.shtml

; usage example:
; opt                   = True
; opt@NumberlegendItems = 5
; opt@LegendLabels      = (/"CWRF/ERI", "CWRF", "WRF", "WRFG", "CMM5"/)
; annoid                = create_legend_from_plot(wks_pdf,plot_pdf(0),opt)
; draw(plot_pdf(0))
; frame(wks_pdf)

begin
    ; Make sure we have an XY plot!
    classname = NhlClassName(plot)
    if(classname .ne. "xyPlotClass") then
        print("create_legend_from_plot: input must be an XY plot")
        return(new(1,graphic))
    end if

    lgres = True

    ; Retrieve some resource values from the XY plot.
    getvalues plot@dataspec
        "xyMonoDashPattern"     : mono_dashindex
        "xyMonoMarkLineMode"    : mono_itemtype
        "xyMonoLineColor"       : mono_linecolor
        "xyMonoLineThickness"   : mono_linethickness
        "xyMonoMarkerColor"     : mono_markercolor
        "xyMonoMarker"          : mono_markerindex
        "xyMonoMarkerSize"      : mono_markersize
        "xyMonoMarkerThickness" : mono_markerthickness
    end getvalues

    ; Start constructing the legend resources based on XY plot resource settings.
    getvalues plot@dataspec
        "xyDashPattern"       : DashIndex
        "xyDashPatterns"      : DashIndexes
        "xyMarkLineMode"      : ItemType
        "xyMarkLineModes"     : ItemTypes
        "xyLineColor"         : LineColor
        "xyLineColors"        : LineColors
        "xyLineThicknessF"    : LineThicknessF
        "xyLineThicknesses"   : LineThicknesses
        "xyMarkerColor"       : MarkerColor
        "xyMarkerColors"      : MarkerColors
        "xyMarker"            : MarkerIndex
        "xyMarkers"           : MarkerIndexes
        "xyMarkerSizeF"       : MarkerSizeF
        "xyMarkerSizes"       : MarkerSizes
        "xyMarkerThicknessF"  : MarkerThicknessF
        "xyMarkerThicknesses" : MarkerThicknesses
    end getvalues

    if(mono_dashindex .or. all(ismissing(DashIndexes))) then
        lgres@lgDashIndex     = DashIndex
        lgres@lgMonoDashIndex = True
    else
        lgres@lgDashIndexes   = DashIndexes
        lgres@lgMonoDashIndex = False
    end if
    if(mono_itemtype .or. all(ismissing(ItemTypes))) then
        lgres@lgItemType     = ItemType
        lgres@lgMonoItemType = True
    else
        lgres@lgItemTypes    = ItemTypes
        lgres@lgMonoItemType = False
    end if
    if(mono_linecolor .or. all(ismissing(LineColors))) then
        lgres@lgLineColor     = LineColor
        lgres@lgMonoLineColor = True
    else
        lgres@lgLineColors    = LineColors
        lgres@lgMonoLineColor = False
    end if
    if(mono_linethickness .or. all(ismissing(LineThicknesses))) then
        lgres@lgLineThicknessF    = LineThicknessF
        lgres@lgMonoLineThickness = True
    else
        lgres@lgLineThicknesses   = LineThicknesses
        lgres@lgMonoLineThickness = False
    end if
    if(mono_markercolor .or. all(ismissing(MarkerColors))) then
        lgres@lgMarkerColor     = MarkerColor
        lgres@lgMonoMarkerColor = True
    else
        lgres@lgMarkerColors    = MarkerColors
        lgres@lgMonoMarkerColor = False
    end if
    if(mono_markerindex .or. all(ismissing(MarkerIndexes))) then
        lgres@lgMarkerIndex     = MarkerIndex
        lgres@lgMonoMarkerIndex = True
    else
        lgres@lgMarkerIndexes   = MarkerIndexes
        lgres@lgMonoMarkerIndex = False
    end if
    if(mono_markersize .or. all(ismissing(MarkerSizes))) then
        lgres@lgMarkerSizeF    = MarkerSizeF
        lgres@lgMonoMarkerSize = True
    else
        lgres@lgMarkerSizes    = MarkerSizes
        lgres@lgMonoMarkerSize = False
    end if
    if(mono_markerthickness .or. all(ismissing(MarkerThicknesses))) then
        lgres@lgMarkerThicknessF    = MarkerThicknessF
        lgres@lgMonoMarkerThickness = True
    else
        lgres@lgMarkerThicknesses   = MarkerThicknesses
        lgres@lgMonoMarkerThickness = False
    end if

    ; Get possible labels for legend
    getvalues plot@dataspec
        "xyExplicitLegendLabels" : exp_legend_labels
        "xyExplicitLabels"       : exp_labels
    end getvalues

    ; Figure out which legend labels to use.
    if(.not. all(ismissing(exp_legend_labels))) then
        legend_labels = exp_legend_labels
    else if(.not. all(ismissing(exp_labels))) then
        legend_labels = exp_labels
    else if(opt .and. isatt(opt,"LegendLabels")) then
        legend_labels = opt@LegendLabels
    else
        print("create_legend_from_plot: warning: no labels provided for legend")
        legend_labels = ""
    end if
    end if
    end if

    ; Set the width and height of legend based on XY plot size
    getvalues plot
        "vpWidthF"  : vpw
        "vpHeightF" : vph
    end getvalues
    ; lgres@vpWidthF  = vpw/4.
    ; lgres@vpHeightF = vph/5.
    lgres@vpWidthF  = vpw/3.5
    lgres@vpHeightF = vph/2.8

    lgres@lgPerimOn   = False
    ; lgres@lgItemOrder = (/3, 2, 1, 0, 7, 6, 5, 4/)
    lgres@lgLineThicknessF = 2.0

    ; lgres@lgLabelFontHeightF = 0.05
    lgres@lgLabelFontHeightF   = 0.15  ; 0.1

    ; Get desired number of legend items
    nitems = dimsizes(legend_labels)
    if(opt .and. isatt(opt,"NumberLegendItems")) then
        if(opt@NumberLegendItems .gt. dimsizes(legend_labels)) then
            print("create_legend_from_plot: you don't have enough legend labels specified. Resetting.")
        else
            nitems = opt@NumberLegendItems
        end if
    end if

    ; Create the legend
    lbid = gsn_create_legend(wks, nitems, legend_labels, lgres)  ; create legend

    ; Attach lower right corner of legend to lower right corner of XY plot.
    amres                  = True
    amres@amJust           = "BottomRight"
    ; amres@amParallelPosF   = 0.50
    ; amres@amOrthogonalPosF = 0.50
    amres@amParallelPosF   = -0.20   ; bigger = more right
    amres@amOrthogonalPosF = 0.45    ; bigger = more down

    ; Add the legend to the plot
    annoid = gsn_add_annotation(plot, lbid, amres)  ; add legend to plot

    ; Return the id so it "lives" outside this function.
    return(annoid)
end  ; function create_legend_from_plot


function compute_rainy_day_1d(pr:float, n_day[1]:integer)
; pr(n_day)
begin
    n_rainy_day = 0

    do i_day = 0, n_day - 1
        if (ismissing(pr(i_day)) .or. pr(i_day) .gt. 1.0e5) then
            print("WARNING: missing value encountered, function compute_rainy_day_1d!")
        else
            if (pr(i_day) .gt. 1.0) then
                n_rainy_day = n_rainy_day + 1
            end if
        end if
    end do

    return n_rainy_day
end  ; function compute_rainy_day_1d


function compute_rainy_day_us_cwrf(pr:float, n_day[1]:integer)
; pr(n_day, 138, 195)
; return n_rainy_day(138, 195)
begin
    n_rainy_day = new((/138, 195/), "integer")
    n_rainy_day = 0

    rainy_flag = new((/138, 195/), "integer")

    pr_tmp = new((/n_day, 138, 195/), "float", pr@_FillValue)
    pr_tmp = pr

    if (any(ismissing(pr_tmp))) then
        print("WARNING: missing value involved, function compute_rainy_day_us_cwrf!")
        pr_tmp = where(ismissing(pr_tmp), 0.0, pr_tmp)
    end if

    do i_day = 0, n_day - 1
        rainy_flag = where(pr_tmp(i_day, :, :) .gt. 1.0, 1, 0)
        n_rainy_day = n_rainy_day + rainy_flag
    end do

    return n_rainy_day
end  ; function compute_rainy_day_us_cwrf


function compute_max_dry_spell_1d(pr:float, n_day[1]:integer)
; pr(n_day)
begin
    if (any(ismissing(pr))) then
        print("WARNING: missing value involved, function compute_max_dry_spell_1d!")
    end if

    max_dry_spell_past = 0
    max_dry_spell_curr = 0
    curr_in_dry        = True

    do i_day = 0, n_day - 1
        if (ismissing(pr(i_day))) then
            continue
        end if

        if (pr(i_day) .le. 1.0) then  ; dry
            if (curr_in_dry) then
                max_dry_spell_curr = max_dry_spell_curr + 1
            else
                curr_in_dry        = True
                max_dry_spell_curr = 1
            end if
        else  ; rainy
            if (curr_in_dry) then
                curr_in_dry = False
                if (max_dry_spell_curr .gt. max_dry_spell_past) then
                    max_dry_spell_past = max_dry_spell_curr
                end if
                max_dry_spell_curr = 0
            end if
        end if
    end do 

    ; last dry spell
    if (curr_in_dry .and. (max_dry_spell_curr .gt. max_dry_spell_past)) then
        max_dry_spell_past = max_dry_spell_curr
    end if

    return max_dry_spell_past
end  ; function compute_max_dry_spell_1d


function compute_max_dry_spell_us_cwrf(pr:float, n_day[1]:integer)
; pr(n_day, 138, 195)
; ruturn max_dry_spell(138, 195)
begin
    if (any(ismissing(pr))) then
        print("WARNING: missing value involved, function compute_max_dry_spell_us_cwrf!")
    end if   

    NX = 195
    NY = 138

    max_dry_spell = new((/NY, NX/), "integer")
    max_dry_spell = 0

    do i_y = 0, NY - 1
        do i_x = 0, NX - 1
            max_dry_spell(i_y, i_x) = compute_max_dry_spell_1d(pr(:, i_y, i_x), n_day)
        end do
    end do

    return max_dry_spell
end  ; function compute_max_dry_spell_us_cwrf


function compute_daily_pr_pctl_1d(pr:float, n_day[1]:integer, pctl[1]:float)
; pr(n_day)
begin
    PCTL_RAINY_DAYS_ONLY = False

    if (any(ismissing(pr))) then
        print("WARNING: missing value involved, function compute_daily_pr_pctl_1d!")
    end if

    if (pctl .lt. 0.0 .or. pctl .gt. 1.0) then
        print("ERROR: invalid pctl parameter, function compute_daily_pr_pctl_1d!")
        exit()
    end if

    pr_tmp = new((/n_day/), "float", pr@_FillValue)
    pr_tmp = pr

    if (pctl .eq. 1.0) then
        do i_day = 0, n_day - 1
            if (ismissing(pr_tmp(i_day)) .or. pr_tmp(i_day) .gt. 1.0e5) then
                pr_tmp(i_day) = pr_tmp@_FillValue
            end if
        end do
        pr_pctl = max(pr_tmp)
        
        return pr_pctl
    end if


    do i_day = 0, n_day - 1
        if (ismissing(pr_tmp(i_day)) .or. pr_tmp(i_day) .gt. 1.0e5) then
            pr_tmp(i_day) = 2.0e5
        end if
    end do

    qsort(pr_tmp)

    do end_day = n_day - 1, 0
        if (ismissing(pr_tmp(i_day)) .or. pr_tmp(i_day) .gt. 1.0e5) then
            continue
        else
            break
        end if       
    end do

    ; for PCTL_ONLY_RAINY_DAYS, not in this function
    if (PCTL_RAINY_DAYS_ONLY) then
        do start_day = 0, end_day - 1
            if (pr_tmp(start_day) .gt. 1.0) then
                break
            end if
        end do
    else
        start_day = 0
    end if

    f_day_pctl = start_day + pctl * (end_day - start_day)
    n_day_pctl = floattoint(f_day_pctl)

    if (end_day .eq. n_day_pctl) then
        print ("ERROR: end_day .eq. n_day_pctl, function compute_daily_pr_pctl_1d!")
        exit()
    end if

    pr_pctl = pr_tmp(n_day_pctl) + (f_day_pctl - n_day_pctl) * (pr_tmp(n_day_pctl+1) - pr_tmp(n_day_pctl))

    return pr_pctl
end  ; function compute_daily_pr_pctl_1d


function compute_daily_pr_pctl_us_cwrf(pr:float, n_day[1]:integer, pctl[1]:float)
; pr(n_day, NY, NX)
; return pr_pctl(NY, NX) 
begin
    if (any(ismissing(pr))) then
        print("WARNING: missing value involved, function compute_daily_pr_pctl_us_cwrf!")
    end if

    NX = 195
    NY = 138

    pr_pctl = new((/NY, NX/), "float", pr@_FillValue)

    do i_y = 0, NY - 1
        do i_x = 0, NX - 1
            pr_pctl(i_y, i_x) = compute_daily_pr_pctl_1d(pr(:, i_y, i_x), n_day, pctl)
        end do
    end do 

    return pr_pctl
end  ; function compute_daily_pr_pctl_us_cwrf


function is_in_array(elem, array, n_elem)
; array(n_elem)
begin
    res = False

    do i_elem = 0, n_elem - 1
        if (.not. ismissing(array(i_elem)) .and. elem .eq. array(i_elem)) then 
            res = True 
            break
        end if
    end do

    return res
end  ; function is_in_array


function compute_pr_ets_us_cwrf(pr_obs:float, pr_mod:float, N_ETS[1]:integer  \
    , lcc[138][195]:float, us_landmask[138][195]:float                        \
    , buffer_zone[1]:logical, ocean[1]:logical, us_land_only[1]:logical)
; pr_obs(NY, NX); pr_mod(NY, NX); n_ets: 30 for seasonal mean
; return ets(n_ets)
; if all of buffer_zone, ocean and us_land_only are false, then it's based on all_land
begin
    NX = 195
    NY = 138

    ets = new((/N_ETS/), "float", 1.0e20)

    n_o    = new((/N_ETS/), "float")
 
    n_f    = new((/N_ETS/), "float")
    n_fo   = new((/N_ETS/), "float")
    n_nfo  = new((/N_ETS/), "float")
    n_fno  = new((/N_ETS/), "float")
    e_n_fo = new((/N_ETS/), "float")

    n_o    = 0.0
    n_all  = 0
    
    n_f    = 0.0
    n_fo   = 0.0
    n_fno  = 0.0
    n_nfo  = 0.0
    e_n_fo = 0.0

    do i_y = 0, NY - 1
        do i_x = 0, NX - 1
            if (.not. buffer_zone) then
                if (i_y .le. 13 .or. i_x .le. 13 .or. NY - 14 .le. i_y .or. NX - 14 .le. i_x) then
                    continue
                end if
            end if

            if (.not. ocean) then
                if (lcc(i_y, i_x) .eq. 16) then
                    continue
                end if
            end if

            if (us_land_only) then
                if (us_landmask(i_y, i_x) .eq. 0) then
                    continue
                end if
            end if

            if (ismissing(pr_obs(i_y, i_x)) .or. pr_obs(i_y, i_x) .gt. 1.0e5 .or. ismissing(pr_mod(i_y, i_x)) .or. pr_mod(i_y, i_x) .gt. 1.0e5) then
                continue
            end if


            n_all = n_all + 1

            do i_ets = 0, N_ETS - 1
                if (pr_obs(i_y, i_x) .gt. 0.5*(i_ets+1)) then
                    n_o(i_ets) = n_o(i_ets) + 1
                end if

                if (pr_mod(i_y, i_x) .gt. 0.5*(i_ets+1)) then
                    n_f(i_ets) = n_f(i_ets) + 1
                end if

                if (pr_obs(i_y, i_x) .gt. 0.5*(i_ets+1) .and. pr_mod(i_y, i_x) .gt. 0.5*(i_ets+1)) then
                    n_fo(i_ets) = n_fo(i_ets) + 1
                end if

                if (pr_obs(i_y, i_x) .le. 0.5*(i_ets+1) .and. pr_mod(i_y, i_x) .gt. 0.5*(i_ets+1)) then
                    n_fno(i_ets) = n_fno(i_ets) + 1
                end if

                if (pr_obs(i_y, i_x) .gt. 0.5*(i_ets+1) .and. pr_mod(i_y, i_x) .le. 0.5*(i_ets+1)) then
                    n_nfo(i_ets) = n_nfo(i_ets) + 1
                end if
            end do  ; end of "do i_ets = 0, N_ETS - 1"

        end do  ; end of "do i_x = 0, NX - 1"
    end do      ; end of "do i_y = 0, NY - 1"

    ; print("n_all = " + n_all + ", n_o = " + n_o + ", n_f = " + n_f + ", n_fo = " + n_fo + ", n_fno = " + n_fno + ", n_nfo = " + n_nfo)

    do i_ets = 0, N_ETS - 1
        e_n_fo(i_ets) = n_f(i_ets) * n_o(i_ets) / n_all

        if (n_fno(i_ets) + n_nfo(i_ets) + n_fo(i_ets) - e_n_fo(i_ets) .gt. 0.0) then
            ets(i_ets) = (n_fo(i_ets) - e_n_fo(i_ets)) / (n_fno(i_ets) + n_nfo(i_ets) + n_fo(i_ets) - e_n_fo(i_ets))
        else
            ets(i_ets) = 0.0
        end if
    end do  ; end of "do i_ets = 0, N_ETS - 1"

    return ets
end  ; function compute_pr_ets_us_cwrf


procedure define_file_var_att_soda_mom5_5day2monthly(fl_nc_monthly_avg:file)
; local 
begin
    ; define file (global) attributes
    global_att             = True
    global_att@description = "SODA-MOM5 reanalysis Monthly averaged from original 5-day average output, data are on MOM5's grid."
    global_att@contact     = "Ligang Chen: lchen2@umd.edu, UMD/AOSC"
    fileattdef(fl_nc_monthly_avg, global_att)

    ; define file dimensions
    dim_names = (/"time", "xt_ocean", "yt_ocean", "st_ocean", "st_edges_ocean", "xu_ocean", "yu_ocean", "sw_ocean", "sw_edges_ocean"/)
    dim_sizes = (/1     , 1440      , 1070      , 50        , 51              , 1440      , 1070      , 50        , 51              /)  
    dim_unlim = (/True  , False     , False     , False     , False           , False     , False     , False     , False           /)  
    filedimdef(fl_nc_monthly_avg, dim_names, dim_sizes, dim_unlim)

    ; define dimension and other static vars
    filevardef(fl_nc_monthly_avg, "time"          , "double", (/"time"          /)) 
    filevardef(fl_nc_monthly_avg, "xt_ocean"      , "double", (/"xt_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "yt_ocean"      , "double", (/"yt_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "st_ocean"      , "double", (/"st_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "st_edges_ocean", "double", (/"st_edges_ocean"/))
    filevardef(fl_nc_monthly_avg, "xu_ocean"      , "double", (/"xu_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "yu_ocean"      , "double", (/"yu_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "sw_ocean"      , "double", (/"sw_ocean"      /)) 
    filevardef(fl_nc_monthly_avg, "sw_edges_ocean", "double", (/"sw_edges_ocean"/))

;   filevardef(fl_nc_monthly_avg, "geolon_t", "float", (/"yt_ocean", "xt_ocean"/))
;   filevardef(fl_nc_monthly_avg, "geolat_t", "float", (/"yt_ocean", "xt_ocean"/))

;   delete(fl_nc_monthly_avg->xt_ocean@_FillValue)
;   xt_ocean_att                = 0
;   xt_ocean_att@_FillValue     = 
;   xt_ocean_att@long_name      = "tcell longitude"
;   xt_ocean_att@units          = "degrees_E"
;   xt_ocean_att@cartesian_axis = "X"
;   filevarattdef(fl_nc_monthly_avg, "xt_ocean", xt_ocean_att)

    
    fl_nc_monthly_avg->time = 1 
    time_att                = True 
    time_att@long_name      = "time"
    time_att@units          = "month"
    time_att@cartesian_axis = "T"
    time_att@calendar       = "JULIAN"
    time_att@bounds         = "time_bounds"
    filevarattdef(fl_nc_monthly_avg, "time", time_att)



    ; define file variables with their dimensions and attributes
    filevardef(fl_nc_monthly_avg, "temp"     , "float", (/"time", "st_ocean", "yt_ocean", "xt_ocean"/))
    filevardef(fl_nc_monthly_avg, "salt"     , "float", (/"time", "st_ocean", "yt_ocean", "xt_ocean"/))
    filevardef(fl_nc_monthly_avg, "u"        , "float", (/"time", "st_ocean", "yu_ocean", "xu_ocean"/))
    filevardef(fl_nc_monthly_avg, "v"        , "float", (/"time", "st_ocean", "yu_ocean", "xu_ocean"/))
    filevardef(fl_nc_monthly_avg, "sea_level", "float", (/"time",             "yt_ocean", "xt_ocean"/))
    filevardef(fl_nc_monthly_avg, "wt"       , "float", (/"time", "sw_ocean", "yt_ocean", "xt_ocean"/))
    filevardef(fl_nc_monthly_avg, "tau_x"    , "float", (/"time",             "yu_ocean", "xu_ocean"/))
    filevardef(fl_nc_monthly_avg, "tau_y"    , "float", (/"time",             "yu_ocean", "xu_ocean"/))

    fl_nc_monthly_avg->temp = 0.0

    temp_att               = True
    temp_att@long_name     = "Potential temperature"
    temp_att@units         = "degrees C"
    temp_att@valid_range   = "-10.f, 500.f"  ; may need ""
;   temp_att@valid_range   = -10., 500.      ; may need ""
    temp_att@missing_value = -1.e+20  ; shouldn't have "f"
;   temp_att@_FillValue    = -1.e+20
;   temp_att@_FillValue    = temp_att@missing_value
    temp_att@cell_methods  = "time: mean"
    temp_att@coordinates   = "geolon_t geolat_t"
    temp_att@standard_name = "sea_water_potential_temperature"
    filevarattdef(fl_nc_monthly_avg, "temp", temp_att)

    salt_att               = True
    salt_att@long_name     = "Practical Salinity"
    salt_att@units         = "psu"
    salt_att@valid_range   = "-10.f, 100.f"  ; may need ""
    salt_att@missing_value = -1.e+20
    salt_att@_FillValue    = -1.e+20
    salt_att@cell_methods  = "time: mean"
    salt_att@coordinates   = "geolon_t geolat_t"
    salt_att@standard_name = "sea_water_salinity"
    filevarattdef(fl_nc_monthly_avg, "salt", salt_att)

    u_att               = True
    u_att@long_name     = "i-current"
    u_att@units         = "m/sec"
    u_att@valid_range   = "-10.f, 10.f"  ; may need ""
    u_att@missing_value = -1.e+20
    u_att@_FillValue    = -1.e+20
    u_att@cell_methods  = "time: mean"
    u_att@coordinates   = "geolon_c geolat_c"
    u_att@standard_name = "sea_water_x_velocity"
    filevarattdef(fl_nc_monthly_avg, "u", u_att)

    v_att               = True
    v_att@long_name     = "j-current"
    v_att@units         = "m/sec"
    v_att@valid_range   = "-10.f, 10.f"  ; may need ""
    v_att@missing_value = -1.e+20
    v_att@_FillValue    = -1.e+20
    v_att@cell_methods  = "time: mean"
    v_att@coordinates   = "geolon_c geolat_c"
    v_att@standard_name = "sea_water_y_velocity"
    filevarattdef(fl_nc_monthly_avg, "v", v_att)

    sea_level_att               = True
    sea_level_att@long_name     = "effective sea level (eta_t + patm/(rho0*g)) on T cells"
    sea_level_att@units         = "meter"
    sea_level_att@valid_range   = "-1000.f, 1000.f"  ; may need ""
    sea_level_att@missing_value = -1.e+20
    sea_level_att@_FillValue    = -1.e+20
    sea_level_att@cell_methods  = "time: mean"
    sea_level_att@coordinates   = "geolon_t geolat_t"
    sea_level_att@standard_name = "sea_surface_height_above_geoid"
    filevarattdef(fl_nc_monthly_avg, "sea_level", sea_level_att)

    wt_att               = True
    wt_att@long_name     = "dia-surface velocity T-points"
    wt_att@units         = "m/sec"
    wt_att@valid_range   = "-100000.f, 100000.f"  ; may need ""
    wt_att@missing_value = -1.e+20
    wt_att@_FillValue    = -1.e+20
    wt_att@cell_methods  = "time: mean"
    wt_att@coordinates   = "geolon_t geolat_t"
    filevarattdef(fl_nc_monthly_avg, "wt", wt_att)

    tau_x_att               = True
    tau_x_att@long_name     = "i-directed wind stress forcing u-velocity"
    tau_x_att@units         = "N/m^2"
    tau_x_att@valid_range   = "-10.f, 10.f"
    tau_x_att@missing_value = -1.e+20
    tau_x_att@_FillValue    = -1.e+20
    tau_x_att@cell_methods  = "time: mean"
    tau_x_att@coordinates   = "geolon_c geolat_c"
    tau_x_att@standard_name = "surface_downward_x_stress"
    filevarattdef(fl_nc_monthly_avg, "tau_x", tau_x_att)

    tau_y_att               = True
    tau_y_att@long_name     = "j-directed wind stress forcing v-velocity"
    tau_y_att@units         = "N/m^2"
    tau_y_att@valid_range   = "-10.f, 10.f"
    tau_y_att@missing_value = -1.e+20
    tau_y_att@_FillValue    = -1.e+20
    tau_y_att@cell_methods  = "time: mean"
    tau_y_att@coordinates   = "geolon_c geolat_c"
    tau_y_att@standard_name = "surface_downward_y_stress"
    filevarattdef(fl_nc_monthly_avg, "tau_y", tau_y_att)
end  ; procedure define_file_var_att_soda_mom5_5day2monthly


procedure copy_dim_var_soda_mom5_5day2monthly(fl_nc_monthly_avg:file, fl_nc_dim:file)
; local
begin
    fl_nc_monthly_avg->xt_ocean       = fl_nc_dim->xt_ocean
    fl_nc_monthly_avg->yt_ocean       = fl_nc_dim->yt_ocean
    fl_nc_monthly_avg->st_ocean       = fl_nc_dim->st_ocean
    fl_nc_monthly_avg->st_edges_ocean = fl_nc_dim->st_edges_ocean
    fl_nc_monthly_avg->xu_ocean       = fl_nc_dim->xu_ocean
    fl_nc_monthly_avg->yu_ocean       = fl_nc_dim->yu_ocean
    fl_nc_monthly_avg->sw_ocean       = fl_nc_dim->sw_ocean
    fl_nc_monthly_avg->sw_edges_ocean = fl_nc_dim->sw_edges_ocean
    fl_nc_monthly_avg->geolon_t       = fl_nc_dim->geolon_t
    fl_nc_monthly_avg->geolat_t       = fl_nc_dim->geolat_t   
end  ; procedure copy_dim_var_soda_mom5_5day2monthly(fl_nc_monthly_avg:file, fl_nc_dim:file)


procedure assign_monthly_avg_soda_mom5_5day2monthly(fl_nc_monthly_avg:file, avg_temp[1][50][1070][1440]:float  \
    , avg_salt[1][50][1070][1440]:float \
    , avg_u:float, avg_v:float, avg_sea_level:float, avg_wt:float, avg_tau_x:float, avg_tau_x:float)
; local
begin
    str_var_names = getfilevarnames(fl_nc_monthly_avg)
    print(str_var_names)

;   str_temp_att = getfilevaratts (fl_nc_monthly_avg, "temp")
;   str_temp_typ = getfilevartypes(fl_nc_monthly_avg, "temp")
;   temp_fv      = getFillValue   (fl_nc_monthly_avg->temp)
;   avg_temp_fv  = getFillValue(avg_temp)
;   print(str_temp_att)
;   print(str_temp_typ)
;   print(temp_fv )
;   print(avg_temp_fv)


    fl_nc_monthly_avg->temp(:, :, :, :)      = avg_temp(:, :, :, :)
    fl_nc_monthly_avg->salt      = avg_salt
    fl_nc_monthly_avg->u         = avg_u
    fl_nc_monthly_avg->v         = avg_v
    fl_nc_monthly_avg->sea_level = avg_sea_level
    fl_nc_monthly_avg->wt        = avg_wt
    fl_nc_monthly_avg->tau_x     = avg_tau_x
    fl_nc_monthly_avg->tau_y     = avg_tau_y
end  ;  procedure assign_monthly_avg_soda_mom5_5day2monthly


; compute the julian day corresponding to the day on the gregorian calendar from SODA
function soda_jday(mon, day, year)
begin
    res_jday = 0




    return res_index
end  ; function get_record_index
