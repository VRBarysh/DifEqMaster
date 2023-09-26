object Dump25DEl: TDump25DEl
  Left = 192
  Top = 106
  Width = 800
  Height = 600
  Caption = 'Dump 25D Electron model'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  OnPaint = FormPaint
  PixelsPerInch = 96
  TextHeight = 13
  object TrackBar1: TTrackBar
    Left = 0
    Top = 548
    Width = 737
    Height = 17
    Max = 10000
    Orientation = trHorizontal
    Frequency = 1
    Position = 0
    SelEnd = 0
    SelStart = 0
    TabOrder = 0
    ThumbLength = 10
    TickMarks = tmBottomRight
    TickStyle = tsNone
    OnChange = TrackBar1Change
  end
  object Edit_t: TLabeledEdit
    Left = 664
    Top = 16
    Width = 121
    Height = 21
    EditLabel.Width = 27
    EditLabel.Height = 13
    EditLabel.Caption = 'Edit_t'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 1
    Text = '0'
  end
  object Edit_nXSteps: TLabeledEdit
    Left = 663
    Top = 56
    Width = 121
    Height = 21
    EditLabel.Width = 40
    EditLabel.Height = 13
    EditLabel.Caption = 'nXSteps'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 2
    Text = '100'
  end
  object Edit_dx: TLabeledEdit
    Left = 663
    Top = 96
    Width = 121
    Height = 21
    EditLabel.Width = 11
    EditLabel.Height = 13
    EditLabel.Caption = 'dx'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 3
    Text = '0.1'
  end
  object Button_XZSlice: TButton
    Left = 704
    Top = 520
    Width = 75
    Height = 25
    Caption = 'XZSlice'
    TabOrder = 4
    OnClick = Button_XZSliceClick
  end
  object Button_XTSlice: TButton
    Left = 704
    Top = 488
    Width = 75
    Height = 25
    Caption = 'XTSlice'
    TabOrder = 5
    OnClick = Button_XTSliceClick
  end
  object Edit_pos1: TLabeledEdit
    Left = 664
    Top = 136
    Width = 57
    Height = 21
    EditLabel.Width = 24
    EditLabel.Height = 13
    EditLabel.Caption = 'Pos1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 6
    Text = '250'
  end
  object Edit_scale1: TLabeledEdit
    Left = 728
    Top = 136
    Width = 57
    Height = 21
    EditLabel.Width = 33
    EditLabel.Height = 13
    EditLabel.Caption = 'Scale1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 7
    Text = '35'
  end
  object Edit_pos2: TLabeledEdit
    Left = 664
    Top = 176
    Width = 57
    Height = 21
    EditLabel.Width = 24
    EditLabel.Height = 13
    EditLabel.Caption = 'Pos2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 8
    Text = '450'
  end
  object Edit_scale2: TLabeledEdit
    Left = 728
    Top = 176
    Width = 57
    Height = 21
    EditLabel.Width = 33
    EditLabel.Height = 13
    EditLabel.Caption = 'Scale2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 9
    Text = '35'
  end
  object Edit_scale3: TLabeledEdit
    Left = 728
    Top = 216
    Width = 57
    Height = 21
    EditLabel.Width = 33
    EditLabel.Height = 13
    EditLabel.Caption = 'Scale3'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 10
    Text = '100'
  end
  object Edit_angle1: TLabeledEdit
    Left = 664
    Top = 256
    Width = 57
    Height = 21
    EditLabel.Width = 33
    EditLabel.Height = 13
    EditLabel.Caption = 'Angle1'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 11
    Text = '0.5'
  end
  object Edit_angle2: TLabeledEdit
    Left = 728
    Top = 256
    Width = 57
    Height = 21
    EditLabel.Width = 33
    EditLabel.Height = 13
    EditLabel.Caption = 'Angle2'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 12
    Text = '0.5'
  end
  object Button_ZTSlice: TButton
    Left = 704
    Top = 456
    Width = 75
    Height = 25
    Caption = 'ZTSlice'
    TabOrder = 13
    OnClick = Button_ZTSliceClick
  end
  object Check_DoubleX: TCheckBox
    Left = 664
    Top = 288
    Width = 97
    Height = 17
    Caption = 'Double X planes'
    TabOrder = 14
  end
  object Check_DrawLines: TCheckBox
    Left = 664
    Top = 312
    Width = 97
    Height = 17
    Caption = 'Draw Lines'
    TabOrder = 15
  end
  object Edit_3DSkip: TLabeledEdit
    Left = 663
    Top = 216
    Width = 57
    Height = 21
    EditLabel.Width = 38
    EditLabel.Height = 13
    EditLabel.Caption = '3D Skip'
    LabelPosition = lpAbove
    LabelSpacing = 3
    TabOrder = 16
    Text = '1'
  end
  object Button_Save: TButton
    Left = 704
    Top = 392
    Width = 75
    Height = 25
    Caption = 'Save data'
    TabOrder = 17
    OnClick = Button_SaveClick
  end
  object Button_ZSlice: TButton
    Left = 704
    Top = 424
    Width = 75
    Height = 25
    Caption = 'ZSlice'
    TabOrder = 18
    OnClick = Button_ZSliceClick
  end
end
