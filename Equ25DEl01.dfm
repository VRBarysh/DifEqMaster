object FormEqu25DEl01: TFormEqu25DEl01
  Left = 291
  Top = 67
  Width = 575
  Height = 604
  Caption = '0'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Edit_nSteps: TLabeledEdit
    Left = 0
    Top = 16
    Width = 81
    Height = 21
    EditLabel.Width = 79
    EditLabel.Height = 13
    EditLabel.Caption = 'Number of Steps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 0
    Text = '100'
  end
  object Edit_nEls: TLabeledEdit
    Left = 88
    Top = 16
    Width = 89
    Height = 21
    EditLabel.Width = 96
    EditLabel.Height = 13
    EditLabel.Caption = 'Number of Electrons'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 1
    Text = '32'
  end
  object Edit_tMax: TLabeledEdit
    Left = 88
    Top = 56
    Width = 89
    Height = 21
    EditLabel.Width = 23
    EditLabel.Height = 13
    EditLabel.Caption = 'tMax'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 2
    Text = '10'
  end
  object Edit_L: TLabeledEdit
    Left = 0
    Top = 56
    Width = 81
    Height = 21
    EditLabel.Width = 6
    EditLabel.Height = 13
    EditLabel.Caption = 'L'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 3
    Text = '10'
  end
  object Edit_r: TLabeledEdit
    Left = 0
    Top = 96
    Width = 81
    Height = 21
    EditLabel.Width = 3
    EditLabel.Height = 13
    EditLabel.Caption = 'r'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 4
    Text = '0.01'
  end
  object Edit_rTime: TLabeledEdit
    Left = 88
    Top = 96
    Width = 89
    Height = 21
    EditLabel.Width = 26
    EditLabel.Height = 13
    EditLabel.Caption = 'rTime'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 5
    Text = '1000'
  end
  object Edit_ElTime: TLabeledEdit
    Left = 0
    Top = 136
    Width = 81
    Height = 21
    EditLabel.Width = 35
    EditLabel.Height = 13
    EditLabel.Caption = 'El Time'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 6
    Text = '10'
  end
  object Edit_onSteps: TLabeledEdit
    Left = 0
    Top = 176
    Width = 81
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_onSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 7
    Text = '10'
  end
  object Edit_oV1: TLabeledEdit
    Left = 0
    Top = 216
    Width = 81
    Height = 21
    EditLabel.Width = 43
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_oV1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 8
    Text = '0.1'
  end
  object Edit_oV2: TLabeledEdit
    Left = 88
    Top = 216
    Width = 89
    Height = 21
    EditLabel.Width = 43
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_oV2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 9
    Text = '10'
  end
  object Check_optim: TCheckBox
    Left = 88
    Top = 176
    Width = 89
    Height = 17
    Caption = 'Check_optim'
    TabOrder = 10
  end
  object Check_ElsOnStart: TCheckBox
    Left = 88
    Top = 136
    Width = 89
    Height = 17
    Caption = 'ElsOnStart'
    TabOrder = 11
  end
  object Edit_ICoef: TLabeledEdit
    Left = 0
    Top = 256
    Width = 81
    Height = 21
    EditLabel.Width = 28
    EditLabel.Height = 13
    EditLabel.Caption = 'I Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 12
    Text = '1'
  end
  object Edit_deltaI: TLabeledEdit
    Left = 0
    Top = 296
    Width = 81
    Height = 21
    EditLabel.Width = 31
    EditLabel.Height = 13
    EditLabel.Caption = 'Delta I'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 13
    Text = '0'
  end
  object Edit_deltaS: TLabeledEdit
    Left = 88
    Top = 296
    Width = 89
    Height = 21
    EditLabel.Width = 35
    EditLabel.Height = 13
    EditLabel.Caption = 'Delta S'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 14
    Text = '0'
  end
  object Check_2Phase: TCheckBox
    Left = 88
    Top = 336
    Width = 89
    Height = 17
    Caption = 'DoublePhase'
    TabOrder = 15
  end
  object Edit_AI: TLabeledEdit
    Left = 88
    Top = 256
    Width = 89
    Height = 21
    EditLabel.Width = 10
    EditLabel.Height = 13
    EditLabel.Caption = 'AI'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 16
    Text = '1'
  end
  object Check_Dump: TCheckBox
    Left = 0
    Top = 336
    Width = 81
    Height = 17
    Caption = 'Dump'
    TabOrder = 17
  end
  object Check_Simple2Phase: TCheckBox
    Left = 0
    Top = 360
    Width = 169
    Height = 17
    Caption = 'Use Simple 2Phase Model'
    TabOrder = 18
  end
  object Edit_rndCoef: TLabeledEdit
    Left = 0
    Top = 400
    Width = 81
    Height = 21
    EditLabel.Width = 40
    EditLabel.Height = 13
    EditLabel.Caption = 'rnd Coef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 19
    Text = '0'
  end
  object Edit_Ai_kappa: TLabeledEdit
    Left = 88
    Top = 400
    Width = 89
    Height = 21
    EditLabel.Width = 57
    EditLabel.Height = 13
    EditLabel.Caption = 'Ai for kappa'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 20
    Text = '0'
  end
  object Check_UseEA: TCheckBox
    Left = 0
    Top = 432
    Width = 177
    Height = 17
    Caption = 'Use EA system'
    TabOrder = 21
  end
  object Edit_EAStartTime: TLabeledEdit
    Left = 0
    Top = 472
    Width = 81
    Height = 21
    EditLabel.Width = 65
    EditLabel.Height = 13
    EditLabel.Caption = 'EA Start Time'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 22
    Text = '100'
  end
  object Edit_EAnSteps: TLabeledEdit
    Left = 88
    Top = 512
    Width = 89
    Height = 21
    EditLabel.Width = 50
    EditLabel.Height = 13
    EditLabel.Caption = 'EA nSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 23
    Text = '1'
  end
  object Edit_EADelta: TLabeledEdit
    Left = 0
    Top = 512
    Width = 81
    Height = 21
    EditLabel.Width = 40
    EditLabel.Height = 13
    EditLabel.Caption = 'EA delta'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 24
    Text = '0'
  end
  object Edit_EAEndTime: TLabeledEdit
    Left = 88
    Top = 472
    Width = 89
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'EA End Time'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 25
    Text = '100'
  end
  object Edit_EA_z0: TLabeledEdit
    Left = 0
    Top = 552
    Width = 81
    Height = 21
    EditLabel.Width = 61
    EditLabel.Height = 13
    EditLabel.Caption = 'EA deltaZ z0'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 26
    Text = '100'
  end
  object Edit_EAkz: TLabeledEdit
    Left = 88
    Top = 552
    Width = 89
    Height = 21
    EditLabel.Width = 61
    EditLabel.Height = 13
    EditLabel.Caption = 'EA deltaZ kz'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 27
    Text = '0.1'
  end
  object Edit_nSpeedFracs: TLabeledEdit
    Left = 200
    Top = 16
    Width = 81
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'nSpeedFracs'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 28
    Text = '1'
  end
  object Edit_SpeedDistrFactor: TLabeledEdit
    Left = 296
    Top = 16
    Width = 81
    Height = 21
    EditLabel.Width = 82
    EditLabel.Height = 13
    EditLabel.Caption = 'SpeedDistrFactor'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 29
    Text = '0'
  end
  object Check_VertA: TCheckBox
    Left = 192
    Top = 336
    Width = 89
    Height = 17
    Caption = 'Vertical A'
    TabOrder = 30
  end
  object Check_MoveAData: TCheckBox
    Left = 280
    Top = 336
    Width = 97
    Height = 17
    Caption = 'Move A'
    TabOrder = 31
  end
  object Edit_KappaMode: TLabeledEdit
    Left = 200
    Top = 64
    Width = 81
    Height = 21
    EditLabel.Width = 61
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa Mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 32
    Text = '0'
  end
  object Edit_KappaV1: TLabeledEdit
    Left = 296
    Top = 64
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 33
    Text = '0'
  end
  object Edit_KappaV2: TLabeledEdit
    Left = 200
    Top = 104
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 34
    Text = '0'
  end
  object Edit_KappaV3: TLabeledEdit
    Left = 296
    Top = 104
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V3'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 35
    Text = '0'
  end
  object Check_maxA: TCheckBox
    Left = 152
    Top = 360
    Width = 89
    Height = 17
    Caption = 'Check_maxA'
    TabOrder = 36
  end
  object Check_PSO: TCheckBox
    Left = 264
    Top = 360
    Width = 113
    Height = 17
    Caption = 'PSO Optimization'
    TabOrder = 37
  end
  object Check_LOVModel: TCheckBox
    Left = 360
    Top = 336
    Width = 113
    Height = 17
    Caption = 'LOV Model'
    TabOrder = 38
  end
  object Edit_KappaV4: TLabeledEdit
    Left = 200
    Top = 144
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V4'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 39
    Text = '0'
  end
  object Edit_KappaV5: TLabeledEdit
    Left = 296
    Top = 144
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V5'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 40
    Text = '0'
  end
  object Check_LoadKappaData: TCheckBox
    Left = 192
    Top = 384
    Width = 113
    Height = 17
    Caption = 'LoadKappaData'
    TabOrder = 41
  end
  object Edit_KappaV6: TLabeledEdit
    Left = 200
    Top = 184
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V6'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 42
    Text = '0'
  end
  object Edit_KappaV7: TLabeledEdit
    Left = 296
    Top = 184
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V7'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 43
    Text = '0'
  end
  object Edit_KappaV8: TLabeledEdit
    Left = 200
    Top = 224
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V8'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 44
    Text = '0'
  end
  object Edit_KappaV9: TLabeledEdit
    Left = 296
    Top = 224
    Width = 81
    Height = 21
    EditLabel.Width = 47
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V9'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 45
    Text = '0'
  end
  object Edit_KappaV10: TLabeledEdit
    Left = 200
    Top = 264
    Width = 81
    Height = 21
    EditLabel.Width = 53
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V10'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 46
    Text = '0'
  end
  object Edit_KappaV11: TLabeledEdit
    Left = 296
    Top = 264
    Width = 81
    Height = 21
    EditLabel.Width = 53
    EditLabel.Height = 13
    EditLabel.Caption = 'Kappa V11'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 47
    Text = '0'
  end
  object Check_1stPeak: TCheckBox
    Left = 384
    Top = 360
    Width = 97
    Height = 17
    Caption = '1st Peak'
    TabOrder = 48
  end
  object Edit_Nu: TLabeledEdit
    Left = 200
    Top = 304
    Width = 81
    Height = 21
    EditLabel.Width = 14
    EditLabel.Height = 13
    EditLabel.Caption = 'Nu'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 49
    Text = '0'
  end
  object Edit_ElMode: TLabeledEdit
    Left = 296
    Top = 304
    Width = 81
    Height = 21
    EditLabel.Width = 42
    EditLabel.Height = 13
    EditLabel.Caption = 'El. Mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 50
    Text = '0'
  end
  object Edit_DeltaV1: TLabeledEdit
    Left = 200
    Top = 552
    Width = 81
    Height = 21
    EditLabel.Width = 38
    EditLabel.Height = 13
    EditLabel.Caption = 'DeltaV1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 51
    Text = '0'
  end
  object Edit_DeltaV2: TLabeledEdit
    Left = 296
    Top = 552
    Width = 81
    Height = 21
    EditLabel.Width = 38
    EditLabel.Height = 13
    EditLabel.Caption = 'DeltaV2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 52
    Text = '0'
  end
  object Edit_AInputMode: TLabeledEdit
    Left = 200
    Top = 472
    Width = 81
    Height = 21
    EditLabel.Width = 61
    EditLabel.Height = 13
    EditLabel.Caption = 'AInput Mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 53
    Text = '0'
  end
  object Edit_AInputV1: TLabeledEdit
    Left = 296
    Top = 472
    Width = 81
    Height = 21
    EditLabel.Width = 44
    EditLabel.Height = 13
    EditLabel.Caption = 'AInputV1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 54
    Text = '0'
  end
  object Edit_AInputV3: TLabeledEdit
    Left = 296
    Top = 512
    Width = 81
    Height = 21
    EditLabel.Width = 44
    EditLabel.Height = 13
    EditLabel.Caption = 'AInputV3'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 55
    Text = '0'
  end
  object Edit_AInputV2: TLabeledEdit
    Left = 200
    Top = 512
    Width = 81
    Height = 21
    EditLabel.Width = 44
    EditLabel.Height = 13
    EditLabel.Caption = 'AInputV2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 56
    Text = '0'
  end
  object Edit_DeltaZMode: TLabeledEdit
    Left = 384
    Top = 472
    Width = 81
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'DeltaZ Mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 57
    Text = '0'
  end
  object Edit_nDeltaZSteps: TLabeledEdit
    Left = 480
    Top = 472
    Width = 81
    Height = 21
    EditLabel.Width = 65
    EditLabel.Height = 13
    EditLabel.Caption = 'nDeltaZSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 58
    Text = '0'
  end
  object Edit_DeltaZMax: TLabeledEdit
    Left = 384
    Top = 512
    Width = 81
    Height = 21
    EditLabel.Width = 52
    EditLabel.Height = 13
    EditLabel.Caption = 'DeltaZMax'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 59
    Text = '0'
  end
  object LabeledEdit4: TLabeledEdit
    Left = 480
    Top = 512
    Width = 81
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'LabeledEdit1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 60
  end
  object LabeledEdit5: TLabeledEdit
    Left = 384
    Top = 552
    Width = 81
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'LabeledEdit1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 61
  end
  object LabeledEdit6: TLabeledEdit
    Left = 480
    Top = 552
    Width = 81
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'LabeledEdit1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 62
  end
  object Check_LoadDeltaZData: TCheckBox
    Left = 192
    Top = 408
    Width = 113
    Height = 17
    Caption = 'LoadDeltaZData'
    TabOrder = 63
  end
  object Edit_TSliceSave: TLabeledEdit
    Left = 392
    Top = 16
    Width = 81
    Height = 21
    EditLabel.Width = 55
    EditLabel.Height = 13
    EditLabel.Caption = 'TSliceSave'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 64
    Text = '-1'
  end
  object Edit_TSliceSaveZ: TLabeledEdit
    Left = 480
    Top = 16
    Width = 81
    Height = 21
    EditLabel.Width = 62
    EditLabel.Height = 13
    EditLabel.Caption = 'TSliceSaveZ'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 65
    Text = '-1'
  end
  object Edit_nThreads: TLabeledEdit
    Left = 312
    Top = 400
    Width = 89
    Height = 21
    EditLabel.Width = 45
    EditLabel.Height = 13
    EditLabel.Caption = 'nThreads'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 66
    Text = '0'
  end
  object Edit_ACoef: TLabeledEdit
    Left = 480
    Top = 64
    Width = 81
    Height = 21
    EditLabel.Width = 29
    EditLabel.Height = 13
    EditLabel.Caption = 'ACoef'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 67
    Text = '0'
  end
  object Edit_RandSeed: TLabeledEdit
    Left = 480
    Top = 104
    Width = 81
    Height = 21
    EditLabel.Width = 51
    EditLabel.Height = 13
    EditLabel.Caption = 'RandSeed'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 68
    Text = '0'
  end
end
