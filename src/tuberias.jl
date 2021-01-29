module tuberias
    
export tuberia
using Unitful

const NPS5 = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 30]
const S5i = [18, 23.4, 30.1, 38.9, 45, 57, 68.78, 84.68, 97.38, 110.08, 135.76, 162.76, 213.56, 266.2, 315.88, 347.68, 398.02, 448.62, 498.44, 549.44, 598.92, 749.3]u"mm"
const S5o = [21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 762]u"mm"
const S5t = [1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 2.11, 2.11, 2.11, 2.11, 2.77, 2.77, 2.77, 3.4, 3.96, 3.96, 4.19, 4.19, 4.78, 4.78, 5.54, 6.35]u"mm"

const NPS10 = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]
const S10i = [7.82, 10.4, 13.8, 17.08, 22.48, 27.86, 36.66, 42.76, 54.76, 66.9, 82.8, 95.5, 108.2, 134.5, 161.5, 211.58, 264.62, 314.66, 342.9, 393.7, 444.3, 495.3, 546.3, 597.3, 644.16, 695.16, 746.16, 797.16, 848.16, 898.16]u"mm"
const S10o = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 660, 711, 762, 813, 864, 914]u"mm"
const S10t = [1.24, 1.65, 1.65, 2.11, 2.11, 2.77, 2.77, 2.77, 2.77, 3.05, 3.05, 3.05, 3.05, 3.4, 3.4, 3.76, 4.19, 4.57, 6.35, 6.35, 6.35, 6.35, 6.35, 6.35, 7.92, 7.92, 7.92, 7.92, 7.92, 7.92]u"mm"

const NPS20 = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]
const S20i = [206.4, 260.3, 311.1, 339.76, 390.56, 441.16, 488.94, 539.94, 590.94, 634.6, 685.6, 736.6, 787.6, 838.6, 888.6]u"mm"
const S20o = [219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 660, 711, 762, 813, 864, 914]u"mm"
const S20t = [6.35, 6.35, 6.35, 7.92, 7.92, 7.92, 9.53, 9.53, 9.53, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7]u"mm"

const NPS30 = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 30, 32, 34, 36]
const S30i = [7.4, 10, 13.4, 16.48, 21.88, 27.6, 36.26, 41.94, 53.94, 63.44, 79.34, 92.04, 104.74, 205.02, 257.4, 307.04, 336.54, 387.34, 434.74, 482.6, 533.6, 581.46, 679.24, 730.24, 781.24, 832.24, 882.24]u"mm"
const S30o = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 711, 762, 813, 864, 914]u"mm"
const S30t = [1.45, 1.85, 1.85, 2.41, 2.41, 2.9, 2.97, 3.18, 3.18, 4.78, 4.78, 4.78, 4.78, 7.04, 7.8, 8.38, 9.53, 9.53, 11.13, 12.7, 12.7, 14.27, 15.88, 15.88, 15.88, 15.88, 15.88]u"mm"

const NPS40 = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 24, 32, 34, 36]
const S40i = [6.84, 9.22, 12.48, 15.76, 20.96, 26.64, 35.08, 40.94, 52.48, 62.68, 77.92, 90.12, 102.26, 128.2, 154.08, 202.74, 254.46, 303.18, 333.34, 381, 428.46, 477.82, 575.04, 778.04, 829.04, 875.9]u"mm"
const S40o = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 610, 813, 864, 914]u"mm"
const S40t = [1.73, 2.24, 2.31, 2.77, 2.87, 3.38, 3.56, 3.68, 3.91, 5.16, 5.49, 5.74, 6.02, 6.55, 7.11, 8.18, 9.27, 10.31, 11.13, 12.7, 14.27, 15.09, 17.48, 17.48, 17.48, 19.05]u"mm"

const NPS60 = [8, 10, 12, 14, 16, 18, 20, 22, 24]
const S60i = [198.48, 247.6, 295.26, 325.42, 373.08, 418.9, 466.76, 514.54, 560.78]u"mm"
const S60o = [219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S60t = [10.31, 12.7, 14.27, 15.09, 16.66, 19.05, 20.62, 22.23, 24.61]u"mm"

const NPS80 = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
const S80i = [5.48, 7.66, 10.7, 13.84, 18.88, 24.3, 32.5, 38.14, 49.22, 58.98, 73.66, 85.44, 97.18, 122.24, 146.36, 193.7, 242.82, 288.84, 317.5, 363.52, 409.34, 455.62, 501.84, 548.08]u"mm"
const S80o = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S80t = [2.41, 3.02, 3.2, 3.73, 3.91, 4.55, 4.85, 5.08, 5.54, 7.01, 7.62, 8.08, 8.56, 9.53, 10.97, 12.7, 15.09, 17.48, 19.05, 21.44, 23.83, 26.19, 28.58, 30.96]u"mm"

const NPS100 = [8, 10, 12, 14, 16, 18, 20, 22, 24]
const S100i = [188.92, 236.48, 280.92, 307.94, 354.02, 398.28, 442.92, 489.14, 532.22]u"mm"
const S100o = [219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S100t = [15.09, 18.26, 21.44, 23.83, 26.19, 29.36, 32.54, 34.93, 38.89]u"mm"

const NPS120 = [4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
const S120i = [92.04, 115.9, 139.76, 182.58, 230.12, 273, 300.02, 344.48, 387.14, 431.8, 476.44, 517.96]u"mm"
const S120o = [114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S120t = [11.13, 12.7, 14.27, 18.26, 21.44, 25.4, 27.79, 30.96, 34.93, 38.1, 41.28, 46.02]u"mm"

const NPS140 = [8, 10, 12, 14, 16, 18, 20, 22, 24]
const S140i = [177.86, 222.2, 266.64, 292.1, 333.34, 377.66, 419.1, 463.74, 505.26]u"mm"
const S140o = [219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S140t = [20.62, 25.4, 28.58, 31.75, 36.53, 39.67, 44.45, 47.63, 52.37]u"mm"

const NPS160 = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
const S160i = [11.74, 15.58, 20.7, 29.5, 34.02, 42.82, 53.94, 66.64, 87.32, 109.54, 131.78, 173.08, 215.84, 257.16, 284.18, 325.42, 366.52, 407.98, 451.04, 490.92]u"mm"
const S160o = [21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610]u"mm"
const S160t = [4.78, 5.56, 6.35, 6.35, 7.14, 8.74, 9.53, 11.13, 13.49, 15.88, 18.26, 23.01, 28.58, 33.32, 35.71, 40.49, 45.24, 50.01, 53.98, 59.54]u"mm"

# Schedules designated STD, XS, and XXS
const NPSSTD = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48]
const STDi = [6.84, 9.22, 12.48, 15.76, 20.96, 26.64, 35.08, 40.94, 52.48, 62.68, 77.92, 90.12, 102.26, 128.2, 154.08, 202.74, 254.46, 304.74, 336.54, 387.34, 437.94, 488.94, 539.94, 590.94, 640.94, 691.94, 742.94, 793.94, 844.94, 894.94, 945.94, 996.94, 1047.94, 1098.94, 1148.94, 1199.94]u"mm"
const STDo = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 660, 711, 762, 813, 864, 914, 965, 1016, 1067, 1118, 1168, 1219]u"mm"
const STDt = [1.73, 2.24, 2.31, 2.77, 2.87, 3.38, 3.56, 3.68, 3.91, 5.16, 5.49, 5.74, 6.02, 6.55, 7.11, 8.18, 9.27, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53, 9.53]u"mm"

const NPSXS = [0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48]
const XSi = [5.48, 7.66, 10.7, 13.84, 18.88, 24.3, 32.5, 38.14, 49.22, 58.98, 73.66, 85.44, 97.18, 122.24, 146.36, 193.7, 247.6, 298.4, 330.2, 381, 431.6, 482.6, 533.6, 584.6, 634.6, 685.6, 736.6, 787.6, 838.6, 888.6, 939.6, 990.6, 1041.6, 1092.6, 1142.6, 1193.6]u"mm"
const XSo = [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 101.6, 114.3, 141.3, 168.3, 219.1, 273, 323.8, 355.6, 406.4, 457, 508, 559, 610, 660, 711, 762, 813, 864, 914, 965, 1016, 1067, 1118, 1168, 1219]u"mm"
const XSt = [2.41, 3.02, 3.2, 3.73, 3.91, 4.55, 4.85, 5.08, 5.54, 7.01, 7.62, 8.08, 8.56, 9.53, 10.97, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7]u"mm"

const NPSXXS = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12]
const XXSi = [6.36, 11.06, 15.22, 22.8, 28, 38.16, 44.96, 58.42, 80.06, 103.2, 124.4, 174.64, 222.2, 273]u"mm"
const XXSo = [21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73, 88.9, 114.3, 141.3, 168.3, 219.1, 273, 323.8]u"mm"
const XXSt = [7.47, 7.82, 9.09, 9.7, 10.15, 11.07, 14.02, 15.24, 17.12, 19.05, 21.95, 22.23, 25.4, 25.4]u"mm"

const schedule_lookup = Dict(     "5"=> (NPS5, S5i, S5o, S5t),
                            "10"=> (NPS10, S10i, S10o, S10t),
                            "20"=> (NPS20, S20i, S20o, S20t),
                            "30"=> (NPS30, S30i, S30o, S30t),
                            "40"=> (NPS40, S40i, S40o, S40t),
                            "60"=> (NPS60, S60i, S60o, S60t),
                            "80"=> (NPS80, S80i, S80o, S80t),
                            "100"=> (NPS100, S100i, S100o, S100t),
                            "120"=> (NPS120, S120i, S120o, S120t),
                            "140"=> (NPS140, S140i, S140o, S140t),
                            "160"=> (NPS160, S160i, S160o, S160t),
                            "STD"=> (NPSSTD, STDi, STDo, STDt),
                            "XS"=> (NPSXS, XSi, XSo, XSt),
                            "XXS"=> (NPSXXS, XXSi, XXSo, XXSt)
)

function Di_lookup(Di, NPSes, Dis, Dos, ts)
    if Dis[end] < Di
        return error("No hay tubería con diámetro interno igual o mayor a este diámetro")
    else
        for i=1:length(Dis) # Go up ascending list; once larger than specified, return
            
            if Dis[i] >= Di
                _nps, _di, _do, _t = NPSes[i], Dis[i], Dos[i], ts[i]
                return (_nps, _di, _do, _t)
            end
        end
    end
end

function Do_lookup(Do, NPSes, Dis, Dos, ts)
    if Dos[end] < Do
        return error("No hay tubería con diámetro externo igual o mayor a este diámetro")
    else
        for i=1:length(Dos) # Go up ascending list; once larger than specified, return
            
            if Dos[i] >= Do
                _nps, _di, _do, _t = NPSes[i], Dis[i], Dos[i], ts[i]
                return (_nps, _di, _do, _t)
            end
        end
    end
end

function NPS_lookup(NPS, NPSes, Dis, Dos, ts)
    i = findfirst(isequal(NPS),NPSes)
    if i==nothing
        error("No hay esa cedula")
    end
    _nps, _di, _do, _t = NPSes[i], Dis[i], Dos[i], ts[i]
    return (_nps, _di, _do, _t)
end

function tuberia(schedule;Do=nothing, Di=nothing, NPS=nothing)
    NPSes, Dis, Dos, ts = schedule_lookup[schedule]
    if Di != nothing
        Di_lookup(Di, NPSes, Dis, Dos, ts)
    elseif Do != nothing
        Do_lookup(Do, NPSes, Dis, Dos, ts)
    elseif NPS != nothing
        NPS_lookup(NPS, NPSes, Dis, Dos, ts)
    end
end 

end

