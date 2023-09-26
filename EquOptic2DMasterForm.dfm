object EquOptic2DForm: TEquOptic2DForm
  Left = 243
  Top = 155
  Width = 450
  Height = 677
  Caption = '2D Optics parameters'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Edit_NStepsZ: TLabeledEdit
    Left = 0
    Top = 16
    Width = 89
    Height = 21
    EditLabel.Width = 40
    EditLabel.Height = 13
    EditLabel.Caption = 'nStepsZ'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 0
    Text = '50'
  end
  object Edit_NStepsX: TLabeledEdit
    Left = 104
    Top = 16
    Width = 89
    Height = 21
    EditLabel.Width = 40
    EditLabel.Height = 13
    EditLabel.Caption = 'nStepsX'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 1
    Text = '50'
  end
  object Edit_TimeP: TLabeledEdit
    Left = 0
    Top = 56
    Width = 89
    Height = 21
    EditLabel.Width = 30
    EditLabel.Height = 13
    EditLabel.Caption = 'TimeP'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 2
    Text = '1'
  end
  object Edit_TimeR: TLabeledEdit
    Left = 104
    Top = 56
    Width = 89
    Height = 21
    EditLabel.Width = 31
    EditLabel.Height = 13
    EditLabel.Caption = 'TimeR'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 3
    Text = '1'
  end
  object Edit_PolNoise: TLabeledEdit
    Left = 0
    Top = 96
    Width = 89
    Height = 21
    EditLabel.Width = 46
    EditLabel.Height = 13
    EditLabel.Caption = 'Pol. noise'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 4
    Text = '0'
  end
  object Edit_tMax: TLabeledEdit
    Left = 104
    Top = 96
    Width = 89
    Height = 21
    EditLabel.Width = 23
    EditLabel.Height = 13
    EditLabel.Caption = 'tMax'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 5
    Text = '10'
  end
  object Edit_A0: TLabeledEdit
    Left = 0
    Top = 520
    Width = 89
    Height = 21
    EditLabel.Width = 13
    EditLabel.Height = 13
    EditLabel.Caption = 'A0'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 6
    Text = '0'
  end
  object Edit_dt: TLabeledEdit
    Left = 104
    Top = 136
    Width = 89
    Height = 21
    EditLabel.Width = 9
    EditLabel.Height = 13
    EditLabel.Caption = 'dt'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 7
    Text = '0.1'
  end
  object Edit_ACoef: TLabeledEdit
    Left = 0
    Top = 136
    Width = 89
    Height = 21
    EditLabel.Width = 32
    EditLabel.Height = 13
    EditLabel.Caption = 'A Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 8
    Text = '0'
  end
  object Edit_FFTSize: TLabeledEdit
    Left = 200
    Top = 560
    Width = 89
    Height = 21
    EditLabel.Width = 42
    EditLabel.Height = 13
    EditLabel.Caption = 'FFT Size'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 9
    Text = '5'
  end
  object Edit_DumpInt: TLabeledEdit
    Left = 0
    Top = 256
    Width = 89
    Height = 21
    EditLabel.Width = 66
    EditLabel.Height = 13
    EditLabel.Caption = 'Dump Interval'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 10
    Text = '1'
  end
  object Check_DumpReq: TCheckBox
    Left = 0
    Top = 288
    Width = 97
    Height = 17
    Caption = 'Dump required'
    TabOrder = 11
  end
  object Edit_nTSubSteps: TLabeledEdit
    Left = 104
    Top = 176
    Width = 89
    Height = 21
    EditLabel.Width = 59
    EditLabel.Height = 13
    EditLabel.Caption = 'nTSubSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 12
    Text = '1'
  end
  object Check_UseSimpleEC: TCheckBox
    Left = 0
    Top = 312
    Width = 97
    Height = 17
    Caption = 'Use Simple EC'
    TabOrder = 13
  end
  object Edit_AGenCoef: TLabeledEdit
    Left = 0
    Top = 176
    Width = 89
    Height = 21
    EditLabel.Width = 49
    EditLabel.Height = 13
    EditLabel.Caption = 'AGenCoef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 14
    Text = '0'
  end
  object Check_UseFastP: TCheckBox
    Left = 0
    Top = 336
    Width = 97
    Height = 17
    Caption = 'Use Fast P'
    Checked = True
    State = cbChecked
    TabOrder = 15
  end
  object Edit_beta: TLabeledEdit
    Left = 104
    Top = 256
    Width = 89
    Height = 21
    EditLabel.Width = 22
    EditLabel.Height = 13
    EditLabel.Caption = 'Beta'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 16
    Text = '0'
  end
  object Check_LockX: TCheckBox
    Left = 0
    Top = 360
    Width = 97
    Height = 17
    Caption = 'Lock X Borders'
    TabOrder = 17
  end
  object Edit_Q: TLabeledEdit
    Left = 104
    Top = 296
    Width = 89
    Height = 21
    EditLabel.Width = 8
    EditLabel.Height = 13
    EditLabel.Caption = 'Q'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 18
    Text = '0'
  end
  object Edit_Z0RefCoef: TLabeledEdit
    Left = 104
    Top = 336
    Width = 89
    Height = 21
    EditLabel.Width = 52
    EditLabel.Height = 13
    EditLabel.Caption = 'Z Ref Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 19
    Text = '0'
  end
  object Edit_A1: TLabeledEdit
    Left = 104
    Top = 520
    Width = 89
    Height = 21
    EditLabel.Width = 13
    EditLabel.Height = 13
    EditLabel.Caption = 'A1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 20
    Text = '0'
  end
  object Edit_nStepsZB: TLabeledEdit
    Left = 0
    Top = 560
    Width = 89
    Height = 21
    EditLabel.Width = 71
    EditLabel.Height = 13
    EditLabel.Caption = 'nStepsZ Bragg'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 21
    Text = '50'
  end
  object Edit_nStepsZT: TLabeledEdit
    Left = 104
    Top = 560
    Width = 89
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'nSteps Delay'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 22
    Text = '1'
  end
  object Check_Equ01: TCheckBox
    Left = 0
    Top = 384
    Width = 97
    Height = 17
    Caption = 'Equ 01'
    TabOrder = 23
  end
  object Edit_FFTScale: TLabeledEdit
    Left = 296
    Top = 560
    Width = 89
    Height = 21
    EditLabel.Width = 49
    EditLabel.Height = 13
    EditLabel.Caption = 'FFT Scale'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 24
    Text = '0'
  end
  object Edit_ZRefPhase: TLabeledEdit
    Left = 104
    Top = 376
    Width = 89
    Height = 21
    EditLabel.Width = 60
    EditLabel.Height = 13
    EditLabel.Caption = 'Z Ref Phase'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 25
    Text = '0'
  end
  object Edit_Z1RefCoef: TLabeledEdit
    Left = 104
    Top = 416
    Width = 89
    Height = 21
    EditLabel.Width = 58
    EditLabel.Height = 13
    EditLabel.Caption = 'Z1 Ref Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 26
    Text = '0'
  end
  object Edit_DiffractionCoef: TLabeledEdit
    Left = 104
    Top = 456
    Width = 89
    Height = 21
    EditLabel.Width = 48
    EditLabel.Height = 13
    EditLabel.Caption = 'Diff param'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 27
    Text = '0'
  end
  object Check_UseMirrorDelay: TCheckBox
    Left = 0
    Top = 408
    Width = 97
    Height = 17
    Caption = 'Use Mirror Delay'
    TabOrder = 28
  end
  object Edit_A2: TLabeledEdit
    Left = 200
    Top = 520
    Width = 89
    Height = 21
    EditLabel.Width = 13
    EditLabel.Height = 13
    EditLabel.Caption = 'A2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 29
    Text = '0'
  end
  object Edit_nThreads: TLabeledEdit
    Left = 200
    Top = 16
    Width = 89
    Height = 21
    EditLabel.Width = 69
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_nThreads'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 30
    Text = '0'
  end
  object Check_NoMedia: TCheckBox
    Left = 0
    Top = 432
    Width = 97
    Height = 17
    Caption = 'No Media'
    Checked = True
    State = cbChecked
    TabOrder = 31
  end
  object Check_Coaxial: TCheckBox
    Left = 0
    Top = 456
    Width = 97
    Height = 17
    Caption = 'Coaxial'
    TabOrder = 32
  end
  object Edit_MaxE: TLabeledEdit
    Left = 200
    Top = 56
    Width = 89
    Height = 21
    EditLabel.Width = 27
    EditLabel.Height = 13
    EditLabel.Caption = 'MaxE'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 33
    Text = '10000000'
  end
  object Edit_AGenZ: TLabeledEdit
    Left = 0
    Top = 216
    Width = 89
    Height = 21
    EditLabel.Width = 58
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_AGenZ'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 34
    Text = '0'
  end
  object Edit_AGenX: TLabeledEdit
    Left = 104
    Top = 216
    Width = 89
    Height = 21
    EditLabel.Width = 58
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_AGenX'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 35
    Text = '0'
  end
  object Edit_SaveBitmapInterval: TLabeledEdit
    Left = 200
    Top = 96
    Width = 89
    Height = 21
    EditLabel.Width = 92
    EditLabel.Height = 13
    EditLabel.Caption = 'SaveBitmapInterval'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 36
    Text = '0'
  end
  object Edit_BitmapNorm: TLabeledEdit
    Left = 200
    Top = 136
    Width = 89
    Height = 21
    EditLabel.Width = 60
    EditLabel.Height = 13
    EditLabel.Caption = 'Bitmap Norm'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 37
    Text = '0'
  end
  object Edit_StartingInv: TLabeledEdit
    Left = 200
    Top = 296
    Width = 89
    Height = 21
    EditLabel.Width = 54
    EditLabel.Height = 13
    EditLabel.Caption = 'Starting Inv'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 38
    Text = '1'
  end
  object Edit_ABackCoef: TLabeledEdit
    Left = 200
    Top = 216
    Width = 89
    Height = 21
    EditLabel.Width = 78
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_ABackCoef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 39
    Text = '0'
  end
  object Edit_MicrowaveMode: TLabeledEdit
    Left = 200
    Top = 336
    Width = 89
    Height = 21
    EditLabel.Width = 82
    EditLabel.Height = 13
    EditLabel.Caption = 'Microwave Mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 40
    Text = '0'
  end
  object Check_TE: TCheckBox
    Left = 0
    Top = 480
    Width = 105
    Height = 17
    Caption = 'Az only media(TE)'
    TabOrder = 41
  end
  object Edit_PurgeNoiseDelay: TLabeledEdit
    Left = 200
    Top = 176
    Width = 89
    Height = 21
    EditLabel.Width = 88
    EditLabel.Height = 13
    EditLabel.Caption = 'Purge Noise Delay'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 42
    Text = '-1'
  end
  object Check_2wB: TCheckBox
    Left = 112
    Top = 480
    Width = 97
    Height = 17
    Caption = '1D with B'
    TabOrder = 43
  end
  object Edit_Delta: TLabeledEdit
    Left = 200
    Top = 456
    Width = 89
    Height = 21
    EditLabel.Width = 25
    EditLabel.Height = 13
    EditLabel.Caption = 'Delta'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 44
    Text = '0'
  end
  object Check_1DwBMod1: TCheckBox
    Left = 200
    Top = 480
    Width = 89
    Height = 17
    Caption = '1DwB Mod1'
    TabOrder = 45
  end
  object Edit_DeltaCoef: TLabeledEdit
    Left = 200
    Top = 416
    Width = 89
    Height = 21
    EditLabel.Width = 50
    EditLabel.Height = 13
    EditLabel.Caption = 'Delta Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 46
    Text = '0'
  end
  object Edit_nStepsX2: TLabeledEdit
    Left = 296
    Top = 16
    Width = 89
    Height = 21
    EditLabel.Width = 46
    EditLabel.Height = 13
    EditLabel.Caption = 'nStepsX2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 47
    Text = '10'
  end
  object Edit_QGrid: TLabeledEdit
    Left = 200
    Top = 376
    Width = 89
    Height = 21
    EditLabel.Width = 27
    EditLabel.Height = 13
    EditLabel.Caption = 'QGrid'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 48
    Text = '0'
  end
  object Check_ForceLinear: TCheckBox
    Left = 296
    Top = 472
    Width = 97
    Height = 17
    Caption = 'Force Linear'
    TabOrder = 49
  end
  object Edit_ACoef2: TLabeledEdit
    Left = 296
    Top = 56
    Width = 89
    Height = 21
    EditLabel.Width = 38
    EditLabel.Height = 13
    EditLabel.Caption = 'ACoef 2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 50
    Text = '0'
  end
  object Edit_nMediaRows: TLabeledEdit
    Left = 296
    Top = 96
    Width = 89
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'nMediaRows'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 51
    Text = '0'
  end
  object Edit_nMediaRowSize: TLabeledEdit
    Left = 296
    Top = 136
    Width = 89
    Height = 21
    EditLabel.Width = 77
    EditLabel.Height = 13
    EditLabel.Caption = 'nMediaRowSize'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 52
    Text = '8'
  end
  object Edit_APhase: TLabeledEdit
    Left = 296
    Top = 176
    Width = 89
    Height = 21
    EditLabel.Width = 37
    EditLabel.Height = 13
    EditLabel.Caption = 'APhase'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 53
    Text = '0'
  end
  object Edit_nMedia: TLabeledEdit
    Left = 296
    Top = 216
    Width = 89
    Height = 21
    EditLabel.Width = 116
    EditLabel.Height = 13
    EditLabel.Caption = 'nMedia ('#1085#1077#1086#1076#1085#1086#1088#1086#1076#1085#1086#1077')'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 54
    Text = '0'
  end
  object Edit_nMediaGaussCoef: TLabeledEdit
    Left = 296
    Top = 256
    Width = 89
    Height = 21
    EditLabel.Width = 81
    EditLabel.Height = 13
    EditLabel.Caption = 'MediaGaussCoef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 55
    Text = '0'
  end
  object Edit_TimeR2: TLabeledEdit
    Left = 296
    Top = 296
    Width = 89
    Height = 21
    EditLabel.Width = 37
    EditLabel.Height = 13
    EditLabel.Caption = 'TimeR2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 56
    Text = '-1'
  end
  object Check_2braggs: TCheckBox
    Left = 296
    Top = 496
    Width = 121
    Height = 17
    Caption = '2braggs solid media'
    TabOrder = 57
  end
  object Edit_nStepsZB2: TLabeledEdit
    Left = 0
    Top = 600
    Width = 89
    Height = 21
    EditLabel.Width = 80
    EditLabel.Height = 13
    EditLabel.Caption = 'nStepsZ Bragg 2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 58
    Text = '50'
  end
  object Check_3braggs: TCheckBox
    Left = 296
    Top = 520
    Width = 121
    Height = 17
    Caption = '2D1D2D solid media'
    TabOrder = 59
  end
  object Edit_ACoef1D: TLabeledEdit
    Left = 104
    Top = 600
    Width = 89
    Height = 21
    EditLabel.Width = 46
    EditLabel.Height = 13
    EditLabel.Caption = 'ACoef 1D'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 60
    Text = '0'
  end
  object Edit_APhase1D: TLabeledEdit
    Left = 200
    Top = 600
    Width = 89
    Height = 21
    EditLabel.Width = 54
    EditLabel.Height = 13
    EditLabel.Caption = 'APhase 1D'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 61
    Text = '0'
  end
  object Edit_AGenCoef1D: TLabeledEdit
    Left = 296
    Top = 600
    Width = 89
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'AGenCoef1D'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 62
    Text = '0'
  end
end
