object PSOSolverForm: TPSOSolverForm
  Left = 289
  Top = 61
  Width = 290
  Height = 529
  Caption = 'PSOSolverForm'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Edit_deltaReMin: TLabeledEdit
    Left = 8
    Top = 24
    Width = 121
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'Min Re(delta)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 0
    Text = '-1'
  end
  object Edit_deltaReMax: TLabeledEdit
    Left = 152
    Top = 24
    Width = 121
    Height = 21
    EditLabel.Width = 66
    EditLabel.Height = 13
    EditLabel.Caption = 'Max Re(delta)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 1
    Text = '1'
  end
  object Edit_deltaImMin: TLabeledEdit
    Left = 8
    Top = 64
    Width = 121
    Height = 21
    EditLabel.Width = 60
    EditLabel.Height = 13
    EditLabel.Caption = 'Min Im(delta)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 2
    Text = '-1'
  end
  object Edit_deltaImMax: TLabeledEdit
    Left = 152
    Top = 64
    Width = 121
    Height = 21
    EditLabel.Width = 63
    EditLabel.Height = 13
    EditLabel.Caption = 'Max Im(delta)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 3
    Text = '1'
  end
  object Edit_gammaReMin: TLabeledEdit
    Left = 8
    Top = 104
    Width = 121
    Height = 21
    EditLabel.Width = 74
    EditLabel.Height = 13
    EditLabel.Caption = 'Min Re(gamma)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 4
    Text = '-1'
  end
  object Edit_gammaReMax: TLabeledEdit
    Left = 152
    Top = 104
    Width = 121
    Height = 21
    EditLabel.Width = 77
    EditLabel.Height = 13
    EditLabel.Caption = 'Max Re(gamma)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 5
    Text = '1'
  end
  object Edit_gammaImMin: TLabeledEdit
    Left = 8
    Top = 144
    Width = 121
    Height = 21
    EditLabel.Width = 71
    EditLabel.Height = 13
    EditLabel.Caption = 'Min Im(gamma)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 6
    Text = '-1'
  end
  object Edit_gammaImMax: TLabeledEdit
    Left = 152
    Top = 144
    Width = 121
    Height = 21
    EditLabel.Width = 74
    EditLabel.Height = 13
    EditLabel.Caption = 'Max Im(gamma)'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 7
    Text = '1'
  end
  object Edit_lx: TLabeledEdit
    Left = 8
    Top = 192
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'Lx'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 8
    Text = '10'
  end
  object Edit_lz: TLabeledEdit
    Left = 152
    Top = 192
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'Lz'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 9
    Text = '5'
  end
  object Edit_alpha: TLabeledEdit
    Left = 8
    Top = 232
    Width = 121
    Height = 21
    EditLabel.Width = 26
    EditLabel.Height = 13
    EditLabel.Caption = 'alpha'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 10
    Text = '1'
  end
  object Edit_epsilon: TLabeledEdit
    Left = 152
    Top = 232
    Width = 121
    Height = 21
    EditLabel.Width = 34
    EditLabel.Height = 13
    EditLabel.Caption = 'Epsilon'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 11
    Text = '0.1'
  end
  object Edit_nx: TLabeledEdit
    Left = 8
    Top = 272
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'nx'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 12
    Text = '1'
  end
  object Edit_nz: TLabeledEdit
    Left = 152
    Top = 272
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'nz'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 13
    Text = '1'
  end
  object Edit_mode: TLabeledEdit
    Left = 8
    Top = 312
    Width = 121
    Height = 21
    EditLabel.Width = 26
    EditLabel.Height = 13
    EditLabel.Caption = 'mode'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 14
    Text = '0'
  end
  object Edit_LzMin: TLabeledEdit
    Left = 8
    Top = 352
    Width = 121
    Height = 21
    EditLabel.Width = 52
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_LzMin'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 15
    Text = '1'
  end
  object Edit_ScaleSteps: TLabeledEdit
    Left = 152
    Top = 352
    Width = 121
    Height = 21
    EditLabel.Width = 78
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_ScaleSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 16
    Text = '100'
  end
  object Edit_SaveBitmapsInt: TLabeledEdit
    Left = 152
    Top = 312
    Width = 121
    Height = 21
    EditLabel.Width = 98
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_SaveBitmapsInt'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 17
    Text = '0'
  end
end
