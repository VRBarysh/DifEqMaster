object FormOpticBasic: TFormOpticBasic
  Left = 192
  Top = 107
  Width = 696
  Height = 480
  Caption = 'FormOpticBasic'
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
    Left = 8
    Top = 16
    Width = 121
    Height = 21
    EditLabel.Width = 58
    EditLabel.Height = 13
    EditLabel.Caption = 'N of z Steps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 0
    Text = '100'
  end
  object Edit_tMax: TLabeledEdit
    Left = 8
    Top = 56
    Width = 121
    Height = 21
    EditLabel.Width = 26
    EditLabel.Height = 13
    EditLabel.Caption = 'Max t'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 1
    Text = '10'
  end
  object Edit_dt: TLabeledEdit
    Left = 8
    Top = 96
    Width = 121
    Height = 21
    EditLabel.Width = 9
    EditLabel.Height = 13
    EditLabel.Caption = 'dt'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 2
    Text = '0.1'
  end
  object Edit_dz: TLabeledEdit
    Left = 8
    Top = 136
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'dz'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 3
    Text = '0.1'
  end
  object Edit_alpha: TLabeledEdit
    Left = 8
    Top = 176
    Width = 121
    Height = 21
    EditLabel.Width = 26
    EditLabel.Height = 13
    EditLabel.Caption = 'alpha'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 4
    Text = '0.0'
  end
  object Edit_alpha1: TLabeledEdit
    Left = 8
    Top = 216
    Width = 121
    Height = 21
    EditLabel.Width = 32
    EditLabel.Height = 13
    EditLabel.Caption = 'alpha1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 5
    Text = '0.0'
  end
  object Edit_beta: TLabeledEdit
    Left = 8
    Top = 256
    Width = 121
    Height = 21
    EditLabel.Width = 21
    EditLabel.Height = 13
    EditLabel.Caption = 'beta'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 6
    Text = '1.0'
  end
  object Edit_PNoise: TLabeledEdit
    Left = 8
    Top = 296
    Width = 121
    Height = 21
    EditLabel.Width = 82
    EditLabel.Height = 13
    EditLabel.Caption = 'Polarization noise'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 7
    Text = '0.001'
  end
  object CheckBox1: TCheckBox
    Left = 8
    Top = 328
    Width = 97
    Height = 17
    Caption = 'CheckBox1'
    TabOrder = 8
  end
end
