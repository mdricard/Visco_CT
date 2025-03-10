Module Module1
    Public data() As Double, pos() As Double, vel() As Double, tor() As Double, xdata() As Double
    Public emg1(), emg2(), emg3(), fittedArray(), coeffArray(), slope() As Double
    Public TenDegreeRef As Double, nSlopePts, Subject As Integer, cond As String
    Public FileName As String, nChannels As Integer, aCH(15) As Int32
    Public StartDown(30), StartUp(30), TenDeg(30), nDnPts, nUpPts, nCrossings, rep As Integer
    Public EnergyAbsorbed(20), EnergyReturned(20), PeakTorque(20), AveStiffness(20), PeakStiffness(20) As Double
    Public nSeconds As Integer, nSamples As Integer, AnatomicalZero, HorLimbWt, LimbWt As Double
    Public rms1 As Double, rms2 As Double, rms3 As Double
    Public TimeValue, Condition, LastRep As Integer
    Public F1 As Form1  ' Needed to be able to reference objects on Form1
    Public StiffForm As New frmStiff  ' Declare an instance of frmStiff 
    ' Biodex Calibration Wt Torque 53.79792 N-m
    Public Const BiodexPosZero As Double = 0.424909465
    Public Const BiodexTorqueZero As Double = 2.496382141
    Public Const BiodexVelZero As Double = 2.49640096
    Public Const BiodexPosPerVolt As Double = 69.64639961
    Public Const BiodexTorquePerVolt As Double = 297.5475687
    Public Const BiodexVelPerVolt = 208.1396293
    Public Const RadToDegree = 57.295779513

    Public Sub SaveData()
        Dim FileDir As String, timestr As String
        If F1.optPre.Checked = True Then
            timestr = "Pre"
        ElseIf F1.optPost.Checked = True Then
            timestr = "Post"
        ElseIf F1.optPost15.Checked = True Then
            timestr = "Post15"
        ElseIf F1.optPost30.Checked = True Then
            timestr = "Post30"
        ElseIf F1.optPlantarRelaxed.Checked = True Then
            timestr = "Relaxed"
        ElseIf F1.optPlantar10.Checked = True Then
            timestr = "TenNm"
        Else
            MsgBox("You Must Enter Pre or Post Test, Click Save/Read Button After Checking", MsgBoxStyle.Critical, "Check the Pre/Post Condition")
            Exit Sub
        End If
        If F1.optControl.Checked = True Then
            cond = "Control"
        ElseIf F1.optStretch.Checked = True Then
            cond = "Stretch"
        ElseIf F1.optDiathermy.Checked = True Then
            cond = "Diathermy"
        ElseIf F1.optDiaStretch.Checked = True Then
            cond = "DiaStretch"
        ElseIf F1.optTqRelax.Checked = True Then
            cond = "TqRelax"
        ElseIf F1.optFamil.Checked = True Then
            cond = "Familiar"
        Else
            MsgBox("You Must Enter Condition, Click Save/Read Button After Checking", MsgBoxStyle.Critical, "Enter Test Condition")
            Exit Sub
        End If
        FileDir = F1.txtDir.Text
        FileName = FileDir + "S" + F1.txtSubject.Text + cond + timestr + ".txt"
        Dim resp As Integer
        If Dir(FileName) <> "" Then
            resp = MsgBox("Warning File " & FileName & " Already Exists. Do you want to replace it?", MsgBoxStyle.YesNo, "File Warning")
            If resp = vbYes Then
                FileOpen(1, FileName, OpenMode.Output) ' Open file for output.
                Dim i As Integer
                WriteLine(1, nSamples, HorLimbWt)
                For i = 0 To nSamples - 1
                    WriteLine(1, i, emg1(i), emg2(i), emg3(i), tor(i), pos(i), vel(i))
                Next i
                FileClose(1)
            Else
                resp = MsgBox("Change Subject, Trial, Condition Information then click Save Data", MsgBoxStyle.Critical, "Change File Name")
            End If
        Else
            FileOpen(1, FileName, OpenMode.Output) ' Open file for output.
            Dim i As Integer
            WriteLine(1, nSamples, HorLimbWt)
            For i = 0 To nSamples - 1
                WriteLine(1, i, emg1(i), emg2(i), emg3(i), tor(i), pos(i), vel(i))
            Next i
            FileClose(1)
        End If
    End Sub
    Public Sub ReadData()
        Dim FileDir As String, timestr As String
        If F1.optPre.Checked = True Then
            timestr = "Pre"
            TimeValue = 1
        ElseIf F1.optPost.Checked = True Then
            timestr = "Post"
            TimeValue = 2
        ElseIf F1.optPost15.Checked = True Then
            timestr = "Post15"
            TimeValue = 3
        ElseIf F1.optPost30.Checked = True Then
            timestr = "Post30"
            TimeValue = 4
        ElseIf F1.optPlantarRelaxed.Checked = True Then
            timestr = "Relaxed"
            TimeValue = 5
        ElseIf F1.optPlantar10.Checked = True Then
            timestr = "TenNm"
            TimeValue = 6
        Else
            MsgBox("You Must Enter Pre or Post Test, Click Save/Read Button After Checking", MsgBoxStyle.Critical, "Check the Pre/Post Condition")
            Exit Sub
        End If
        If F1.optFamil.Checked = True Then
            cond = "Familiar"
            Condition = 1
        ElseIf F1.optControl.Checked = True Then
            cond = "Control"
            Condition = 2
        ElseIf F1.optStretch.Checked = True Then
            cond = "Stretch"
            Condition = 3
        ElseIf F1.optDiathermy.Checked = True Then
            cond = "Diathermy"
            Condition = 4
        ElseIf F1.optDiaStretch.Checked = True Then
            cond = "DiaStretch"
            Condition = 5
        ElseIf F1.optTqRelax.Checked = True Then
            cond = "TqRelax"
            Condition = 6
        Else
            MsgBox("You Must Enter Condition, Click Save/Read Button After Checking", MsgBoxStyle.Critical, "Enter Test Condition")
            Exit Sub
        End If
        Dim sum As Double
        sum = 0
        FileDir = F1.txtDir.Text
        Try
            FileName = FileDir + "S" + F1.txtSubject.Text + cond + timestr + ".txt"
            FileOpen(1, FileName, OpenMode.Input) ' Open file for input.
            Input(1, nSamples)
            Input(1, HorLimbWt)
            'nSamples = CInt(F1.txtRate.Text) * CInt(F1.txtTime.Text)
            ReDim pos(nSamples - 1), tor(nSamples - 1), vel(nSamples - 1), xdata(nSamples - 1)
            ReDim emg1(nSamples - 1), emg2(nSamples - 1), emg3(nSamples - 1)
            Dim i As Integer, num As Integer
            For i = 0 To nSamples - 1
                Input(1, xdata(i))
                Input(1, emg1(i))
                Input(1, emg2(i))
                Input(1, emg3(i))
                Input(1, tor(i))
                'tor(i) = -1 * tor(i)
                Input(1, pos(i))
                sum += pos(i)
                Input(1, vel(i))
            Next i
            FileClose(1)
            ' Filter Torque, Position & Velocity
            Dim ButterworthLowPassFilter As NationalInstruments.Analysis.Dsp.Filters.ButterworthLowpassFilter = New ButterworthLowpassFilter(4, CDbl(F1.txtRate.Text), CDbl(F1.txtTqFilter.Text))
            tor = ButterworthLowPassFilter.FilterData(tor)
            Dim FilterPosition As NationalInstruments.Analysis.Dsp.Filters.ButterworthLowpassFilter = New ButterworthLowpassFilter(4, CDbl(F1.txtRate.Text), CDbl(F1.txtPsFilter.Text))
            pos = ButterworthLowPassFilter.FilterData(pos)
            vel = ButterworthLowPassFilter.FilterData(vel)

            F1.PlotLeft1.PlotY(tor)
            F1.PlotRight1.PlotY(emg1)
            F1.PlotLeft2.PlotY(pos)
            F1.PlotRight2.PlotY(emg2)
            F1.PlotLeft3.PlotY(vel)
            F1.PlotRight3.PlotY(emg3)
            ComputeRMS(emg1, rms1, nSamples)
            ComputeRMS(emg2, rms2, nSamples)
            ComputeRMS(emg3, rms3, nSamples)
            F1.lblRMS1.Text = CStr(Format(rms1, "#####.#"))
            F1.lblRMS2.Text = CStr(Format(rms2, "#####.#"))
            F1.lblRMS3.Text = CStr(Format(rms3, "#####.#"))

            If F1.chkAnalyze.Checked = True Then
                TenDegreeRef = sum / nSamples
                ZeroCrossings(pos, 100, nSamples - 10, TenDegreeRef, TenDeg, nCrossings)
                LastRep = 12
                For rep = 0 To LastRep - 1
                    GraphLoadAndUnloadCurve(StartUp(rep), StartDown(rep), StartUp(rep + 1))
                Next rep
                '        GraphLoadCurve(StartUp(0), StartDown(0))
                'GraphUnloadCurve(StartUp(0), StartUp(1))
                ' FindEndPoints(pos, 0, nSamples - 1)
                SaveStats()
            End If
        Catch ex As Exception 'When System.IO.FileNotFoundException  ' Catch the error.
            ' MsgBox(ex.ToString)   ' Show friendly error message.'
            MsgBox(FileName + "   Not Found!", MsgBoxStyle.OKOnly, "File Not Found ERROR: Check File Name Settings and Try Again")
        Finally
            Beep()   ' Beep after error processing.
        End Try
    End Sub
    Sub MeanSD(ByVal Curve() As Double, ByVal First As Integer, ByVal Last As Integer, ByRef mean As Double, ByRef sd As Double)
        Dim sum As Double, cnt, i As Integer
        Sum = 0.0#
        cnt = 0
        For i = First To Last
            Sum = Sum + Curve(i)
            cnt = cnt + 1
        Next i
        mean = Sum / cnt
        Sum = 0.0#
        For i = First To Last
            Sum = Sum + (Curve(i) - mean) ^ 2
        Next i
        sd = System.Math.Sqrt(sum / (cnt - 1))
    End Sub

    Public Sub SaveStats()
        Dim sd As Double, i As Integer
        MeanSD(EnergyAbsorbed, 0, LastRep, EnergyAbsorbed(20), sd)
        MeanSD(EnergyReturned, 0, LastRep, EnergyReturned(20), sd)
        MeanSD(PeakTorque, 0, LastRep, PeakTorque(20), sd)
        MeanSD(AveStiffness, 0, LastRep, AveStiffness(20), sd)
        MeanSD(PeakStiffness, 0, LastRep, PeakStiffness(20), sd)
        Dim FileDir As String
        Subject = CInt(F1.txtSubject.Text)
        FileDir = F1.txtDir.Text
        '-------------------------  Save Individual Reps -------------------
        FileName = FileDir + "REPstats.txt"
        If Dir(FileName) = "" Then
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            'EnergyAbsorbed(20), EnergyReturned(20), PeakTorque(20), AveStiffness(20), PeakStiffness(20)
            WriteLine(1, "Subject", "cond", "Time", "rep", "EnergyAbsorbed", "EnergyReturned", "PeakTorque", "AveStiffness", "PeakStiffness")
            For i = 0 To 11
                WriteLine(1, Subject, Condition, TimeValue, i, EnergyAbsorbed(i), EnergyReturned(i), PeakTorque(i), AveStiffness(i), PeakStiffness(i))
            Next i
            FileClose(1)
        Else
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            For i = 0 To 11
                WriteLine(1, Subject, Condition, TimeValue, i, EnergyAbsorbed(i), EnergyReturned(i), PeakTorque(i), AveStiffness(i), PeakStiffness(i))
            Next i
            FileClose(1)
        End If
        '----------------- Save Mean Values ----------------------------
        FileName = FileDir + "AverageStats.txt"
        If Dir(FileName) = "" Then
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            'EnergyAbsorbed(20), EnergyReturned(20), PeakTorque(20), AveStiffness(20), PeakStiffness(20)
            WriteLine(1, "Subject", "cond", "Time", "EnergyAbsorbed", "EnergyReturned", "PeakTorque", "AveStiffness", "PeakStiffness")
            WriteLine(1, Subject, Condition, TimeValue, EnergyAbsorbed(20), EnergyReturned(20), PeakTorque(20), AveStiffness(20), PeakStiffness(20))
            FileClose(1)
        Else
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            WriteLine(1, Subject, Condition, TimeValue, EnergyAbsorbed(20), EnergyReturned(20), PeakTorque(20), AveStiffness(20), PeakStiffness(20))
            FileClose(1)
        End If
        '----------------- Save Mean Values of the Last 3 ----------------------------
        MeanSD(EnergyAbsorbed, LastRep - 2, LastRep, EnergyAbsorbed(19), sd)
        MeanSD(EnergyReturned, LastRep - 2, LastRep, EnergyReturned(19), sd)
        MeanSD(PeakTorque, LastRep - 2, LastRep, PeakTorque(19), sd)
        MeanSD(AveStiffness, LastRep - 2, LastRep, AveStiffness(19), sd)
        MeanSD(PeakStiffness, LastRep - 2, LastRep, PeakStiffness(19), sd)

        FileName = FileDir + "Last3Stats.txt"
        If Dir(FileName) = "" Then
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            'EnergyAbsorbed(19), EnergyReturned(19), PeakTorque(19), AveStiffness(19), PeakStiffness(19)
            WriteLine(1, "Subject", "cond", "Time", "EnergyAbsorbed", "EnergyReturned", "PeakTorque", "AveStiffness", "PeakStiffness")
            WriteLine(1, Subject, Condition, TimeValue, EnergyAbsorbed(19), EnergyReturned(19), PeakTorque(19), AveStiffness(19), PeakStiffness(19))
            FileClose(1)
        Else
            FileOpen(1, FileName, OpenMode.Append) ' Open file for Append.
            WriteLine(1, Subject, Condition, TimeValue, EnergyAbsorbed(19), EnergyReturned(19), PeakTorque(19), AveStiffness(19), PeakStiffness(19))
            FileClose(1)
        End If

    End Sub
    Public Sub ComputeRMS(ByVal Muscle() As Double, ByRef rmsVal As Double, ByVal n As Integer)
        '----------------------------------------------------------------------------------------------------------
        ' Routine computes RMS for an emg signal.
        ' The RMS value is stored in the array rms().
        ' array rms() is dimensioned rms(1 to number of muscles, 1 to number of nTrials)
        ' array EMGdata() is dimensioned EMGdata(1 to number of muscles, 0 to n-1)
        '----------------------------------------------------------------------------------------------------------
        Dim Sum As Double, i As Integer
        Sum = 0.0#
        For i = 0 To n - 1
            Sum = Sum + Muscle(i) ^ 2
        Next i
        rmsVal = 1000 * Sqrt(Sum / 1000.0#)
    End Sub
    Public Sub ZeroCrossings(ByVal Curve() As Double, ByVal FirstPt As Integer, ByVal LastPt As Integer, ByVal RefValue As Double, ByVal Crossings() As Integer, ByRef nCrossings As Integer)
        ' This routine finds the number of times a curve goes from negative to positive
        ' or positive to negative.  The point of zero crossing is the last positive value
        ' or the last negative value, depending upon the direction of sign change.
        ' The variable RefValue is typically set at 0, in this case we will use 10 deg.
        Dim i As Integer
        i = FirstPt
        nUpPts = 0
        nDnPts = 0
        nCrossings = 0
        Do While i < LastPt
            Do While (Curve(i) > RefValue) And (i < LastPt)
                i = i + 1
                If Curve(i) < RefValue Then
                    Crossings(nCrossings) = i - 1
                    ' Now skip 50 pts
                    i = i + 50
                    nCrossings = nCrossings + 1
                End If
            Loop
            Do While (Curve(i) < RefValue) And (i < LastPt)
                i = i + 1
                If Curve(i) > RefValue Then
                    Crossings(nCrossings) = i - 1
                    ' Now skip 50 pts
                    i = i + 50
                    nCrossings = nCrossings + 1
                End If
            Loop
        Loop
        Dim aMax, aMin As Double, amaxPt, aMinPt As Integer
        For i = 0 To nCrossings - 2
            ' Check if the point is going up or down?
            If (pos(Crossings(i) + 50)) < pos(Crossings(i)) Then 'Find the Min
                MaxMin(pos, Crossings(i), Crossings(i + 1), aMax, aMin, amaxPt, StartUp(nUpPts))
                nUpPts += 1
            Else                            'Find the Max
                MaxMin(pos, Crossings(i), Crossings(i + 1), aMax, aMin, StartDown(nDnPts), aMinPt)
                nDnPts += 1
            End If
        Next i
        ' Now Determine the number of repetitions
        'Dim FirstRepStartUpPt, LastRepStartUpPt As Integer
        'i = 0
        'Do Until (StartUp(i) < StartDown(0))
        '    i = i + 1
        '   Loop
        '        FirstRepStartUpPt = i
        '        LastRepStartUpPt = i

    End Sub
    Public Sub GraphLoadAndUnloadCurve(ByVal FirstPt As Integer, ByVal MidPt As Integer, ByVal LastPt As Integer)
        Dim numpoints, i, k As Integer
        Dim Area As Double

        numpoints = MidPt - FirstPt
        Dim OneCurve(numpoints - 1) As Double
        ReDim xdata(numpoints - 1)
        k = 0
        For i = FirstPt To MidPt - 1
            OneCurve(k) = tor(i)
            xdata(k) = pos(i) / RadToDegree
            k += 1
        Next i
        '        StiffForm.Show()
        '        StiffForm.dataPlot.PlotXY(xdata, OneCurve)
        PolyFit(OneCurve, xdata, numpoints, Area, False)   ' This is the area under the load curve

        numpoints = LastPt - MidPt
        Dim xdata2(numpoints - 1), OneCurve2(numpoints - 1) As Double
        k = 0
        For i = MidPt To LastPt - 1
            OneCurve2(k) = tor(i)
            xdata2(k) = pos(i) / RadToDegree
            k += 1
        Next i
        StiffForm.Show()
        'StiffForm.dataPlot.PlotXY(xdata, OneCurve)
        PolyFit(OneCurve2, xdata2, numpoints, EnergyReturned(rep), True)   'Area under the Unload curve
        EnergyAbsorbed(rep) = Area - EnergyReturned(rep)


    End Sub

    Public Sub MaxMin(ByVal Curve() As Double, ByVal First As Integer, ByVal Last As Integer, ByRef Ymax As Double, ByRef Ymin As Double, ByRef MaxT As Integer, ByRef MinT As Integer)
        Dim i As Integer
        Ymax = -10000 : Ymin = 10000
        For i = First To Last - 1
            If Curve(i) > Curve(i + 1) And Curve(i) >= Ymax Then
                Ymax = Curve(i)
                MaxT = i
            ElseIf Curve(i + 1) >= Ymax Then
                Ymax = Curve(i + 1)
                MaxT = i + 1
            End If
        Next i
        For i = First To Last - 1
            If Curve(i) < Curve(i + 1) And Curve(i) <= Ymin Then
                Ymin = Curve(i)
                MinT = i
            ElseIf Curve(i + 1) <= Ymin Then
                Ymin = Curve(i + 1)
                MinT = i + 1
            End If
        Next i
    End Sub
    Public Sub GraphLoadCurve(ByVal FirstPt As Integer, ByVal LastPt As Integer)
        Dim numpoints, i, k As Integer
        Dim Area As Double

        numpoints = LastPt - FirstPt
        Dim OneCurve(numpoints - 1) As Double
        ReDim xdata(numpoints - 1)
        k = 0
        For i = FirstPt To LastPt - 1
            OneCurve(k) = tor(i)
            xdata(k) = pos(i) / RadToDegree
            k += 1
        Next i
        StiffForm.Show()
        StiffForm.dataPlot.PlotXY(xdata, OneCurve)
        PolyFit(OneCurve, xdata, numpoints, False, Area)

    End Sub
    Public Sub GraphUnloadCurve(ByVal FirstPt As Integer, ByVal LastPt As Integer)
        Dim numpoints, i, k As Integer
        Dim Area As Double

        numpoints = LastPt - FirstPt
        Dim OneCurve(numpoints - 1) As Double
        ReDim xdata(numpoints - 1)
        k = 0
        For i = FirstPt To LastPt - 1
            OneCurve(k) = tor(i)
            xdata(k) = pos(i) / RadToDegree
            k += 1
        Next i
        StiffForm.Show()
        StiffForm.dataPlot.PlotXY(xdata, OneCurve)
        PolyFit(OneCurve, xdata, numpoints, False, Area)

        '        StiffForm.XyPointAnnotation1.SetPosition(xdata(StartDown(0) - StartUp(0)), OneCurve(StartDown(0) - StartUp(0)))

    End Sub
    Public Sub PolyFit(ByVal yArray() As Double, ByVal xArray() As Double, ByVal n As Integer, ByRef Area As Double, ByVal IsItUnload As Boolean)
        ' If the Time xArray() has equally spaced points set EqualXIntervals to true
        ' If the Time xArray() does NOT contain EQUALLY SPACED points set EqualXIntervals to false
        ' This Function is Y = 4X^3 from x = 0 to x = 1
        ' The exact area of the function is 1
        ' X^4 evaluated at 0 and 1 yields the area of 1 - 0 = 1
        Dim mean, LowerLimitIntegration, UpperLimitIntegration, IntegralCoeff(4) As Double, i As Integer
        ReDim coeffArray(4)
        fittedArray = NationalInstruments.Analysis.Math.CurveFit.PolynomialFit(xArray, yArray, 4, PolynomialFitAlgorithm.Svd, coeffArray, mean)
        'y = coeffArray(0) + coeffArray(1)*x + coeffArray(2)*x^2 + coeffArray(3)*x^3 + coeffArray(4)*x^4
        'yPrime =  coeffArray(1) + 2*coeffArray(2)*x + 3*coeffArray(3)*x^2 + 4*coeffArray(4)*x^3
        IntegralCoeff(0) = coeffArray(0)
        IntegralCoeff(1) = coeffArray(1) / 2
        IntegralCoeff(2) = coeffArray(2) / 3
        IntegralCoeff(3) = coeffArray(3) / 4
        IntegralCoeff(4) = coeffArray(4) / 5
        LowerLimitIntegration = IntegralCoeff(0) * xArray(0) + (IntegralCoeff(1) * (xArray(0) ^ 2)) + (IntegralCoeff(2) * (xArray(0) ^ 3)) + (IntegralCoeff(3) * (xArray(0) ^ 4)) + (IntegralCoeff(4) * (xArray(0) ^ 5))
        UpperLimitIntegration = IntegralCoeff(0) * xArray(n - 1) + (IntegralCoeff(1) * (xArray(n - 1) ^ 2)) + (IntegralCoeff(2) * (xArray(n - 1) ^ 3)) + (IntegralCoeff(3) * (xArray(n - 1) ^ 4)) + (IntegralCoeff(4) * (xArray(n - 1) ^ 5))
        If IsItUnload = False Then
            Area = UpperLimitIntegration - LowerLimitIntegration
        Else
            ' The unload curve is backwards so the upperlimit is the lowerlimit
            Area = LowerLimitIntegration - UpperLimitIntegration
        End If

        ' StiffForm.FittedPlot.PlotXY(xArray, yArray)
        StiffForm.lblEnergy.Text = Area
        If IsItUnload = False Then    '  Compute slope for Load Curve Only
            ReDim slope(n - 1)
            For i = 0 To n - 1
                'yPrime =  coeffArray(1) + 2*coeffArray(2)*x + 3*coeffArray(3)*x^2 + 4*coeffArray(4)*x^3
                slope(i) = coeffArray(1) + (2 * coeffArray(2) * xArray(i)) + (3 * coeffArray(3) * (xArray(i) ^ 2)) + (4 * coeffArray(4) * (xArray(i) ^ 3))
            Next i
            ' 10 degrees = 0.174533 radians
            Dim xVal, Tor10, ymin As Double, maxt, mint As Integer
            xVal = xArray(n - 1) - 0.174533
            ' Get Torque 10 deg below max angle
            Tor10 = coeffArray(0) + (coeffArray(1) * xVal) + (coeffArray(2) * (xVal ^ 2)) + (coeffArray(3) * (xVal ^ 3)) + (coeffArray(4) * (xVal ^ 4))
            PeakTorque(rep) = coeffArray(0) + (coeffArray(1) * xArray(n - 1)) + (coeffArray(2) * (xArray(n - 1) ^ 2)) + (coeffArray(3) * (xArray(n - 1) ^ 3)) + (coeffArray(4) * (xArray(n - 1) ^ 4))
            AveStiffness(rep) = (PeakTorque(rep) - Tor10) / 0.174533
            MaxMin(slope, 0, n - 1, PeakStiffness(rep), ymin, maxt, mint)

        End If
        'Plot the polynomial function for Load and Unload Curves
        If IsItUnload = False Then    '  Plot the Loading Curve
            StiffForm.Show()
            StiffForm.dataPlot.PlotXY(xArray, fittedArray)
        Else                          '  Plot the UnLoading Curve
            StiffForm.Show()
            StiffForm.FittedPlot.PlotXY(xArray, fittedArray)
            MsgBox("Subject:" + F1.txtSubject.Text, MsgBoxStyle.Critical, "PrePost = " + Str(TimeValue) + "Rep = " + Str(rep))
        End If
    End Sub
    Public Sub OldPolyFit(ByVal yArray() As Double, ByVal xArray() As Double, ByVal n As Integer, ByVal EqualXIntervals As Boolean, ByRef AreaBode10 As Double)
        ' If the Time xArray() has equally spaced points set EqualXIntervals to true
        ' If the Time xArray() does NOT contain EQUALLY SPACED points set EqualXIntervals to false
        ' This Function is Y = 4X^3 from x = 0 to x = 1
        ' The exact area of the function is 1
        ' X^4 evaluated at 0 and 1 yields the area of 1 - 0 = 1
        Dim mean As Double
        Dim SecondDeriv(n - 1), yPrimeFirst, yPrimeLast As Double
        ReDim coeffArray(4)
        fittedArray = NationalInstruments.Analysis.Math.CurveFit.PolynomialFit(xArray, yArray, 4, PolynomialFitAlgorithm.Svd, coeffArray, mean)

        Dim NumberofInterpolationPts As Integer, NewTimeVal As Double
        NumberofInterpolationPts = 100
        Dim yArray10((n * NumberofInterpolationPts) - 1), xArray10((n * NumberofInterpolationPts) - 1) As Double, i As Integer

        Dim Dt, newDt As Double
        If EqualXIntervals = True Then
            Dt = xArray(1) - xArray(0)
            newDt = Dt / NumberofInterpolationPts
        Else
            ' If the X Intervals are not equally spaced then use ave value
            Dt = (xArray(n - 1) - xArray(0)) / n
            newDt = Dt / NumberofInterpolationPts
        End If
        yPrimeFirst = (yArray(1) - yArray(0)) / Dt
        yPrimeLast = (yArray(n - 1) - yArray(n - 2)) / Dt
        SecondDeriv = NationalInstruments.Analysis.Math.CurveFit.SplineInterpolant(xArray, yArray, yPrimeFirst, yPrimeLast)

        NewTimeVal = xArray(0)
        For i = 0 To (n * NumberofInterpolationPts) - 1
            'yArray10(i) = NationalInstruments.Analysis.Math.CurveFit.SplineInterpolation(xArray, yArray, SecondDeriv, NewTimeVal)
            yArray10(i) = NationalInstruments.Analysis.Math.CurveFit.SplineInterpolation(xArray, fittedArray, SecondDeriv, NewTimeVal)
            xArray10(i) = NewTimeVal
            NewTimeVal += newDt
        Next i

        AreaBode10 = NationalInstruments.Analysis.Math.Calculus.NumericIntegration(yArray10, newDt, IntegrationMethod.BodeRule)
        StiffForm.FittedPlot.PlotXY(xArray10, yArray10)
        StiffForm.lblEnergy.Text = AreaBode10
        ' Plot Raw Data and Fitted Data
        ' Plot the data on the graph.
        '   F1.dataPlot.PlotXY(xArray, yArray)
        'F1.fittedPlot.PlotXY(xArray, fittedArray)
        '   F1.fittedPlot.PlotXY(xArray10, yArray10)
        nSlopePts = fittedArray.GetUpperBound(0)

        ReDim slope(nSlopePts - 2)
        For i = 1 To nSlopePts - 2
            slope(i) = (fittedArray(i + 1) - fittedArray(i)) / (xArray(i + 1) - xArray(i))
        Next i

    End Sub
    Public Sub FindEndPoints(ByVal Curve() As Double, ByVal First As Integer, ByVal Last As Integer)
        Dim i As Integer
        nDnPts = 0
        nUpPts = 0
        Dim maxval, minval As Double, minpt, maxpt As Integer
        'nSkipPts = Val(F1.txtSkip.Text)
        i = 100
        Do Until (pos(i) < pos(i + 25)) And (pos(i) < pos(i + 50)) And (pos(i + 50) > pos(i + 75)) And (pos(i) > 25) Or (i = nSamples - 15)
            i = i + 1
        Loop
        MaxMin(pos, i, i + 75, maxval, minval, StartDown(nDnPts), minpt)
        nDnPts = nDnPts + 1
        i = i - 75
        ' Now back up and find the start
        Do Until (pos(i) >= pos(i + 5)) Or (i < 10)
            i = i - 1
        Loop
        StartUp(nUpPts) = i
        nUpPts = nUpPts + 1
        i = StartDown(0) + 100
        Do While i < (nSamples - 12)
            Do Until (pos(i) <= pos(i + 10)) Or (i = (nSamples - 10))
                If i > (nSamples - 20) Then Exit Do
                i = i + 1
            Loop
            StartUp(nUpPts) = i + 2
            nUpPts = nUpPts + 1
            If i > (nSamples - 110) Then
                Exit Do
            Else
                i = i + 100
            End If
            Do Until (pos(i) > pos(i + 10)) Or (i = (nSamples - 10))
                If i > (nSamples - 20) Then Exit Do
                i = i + 1
            Loop
            StartDown(nDnPts) = i + 2
            nDnPts = nDnPts + 1
            i = i + 100
        Loop
    End Sub

End Module
