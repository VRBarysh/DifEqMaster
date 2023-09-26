object Form1: TForm1
  Left = 197
  Top = 122
  Width = 838
  Height = 644
  Caption = 'no'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  OnPaint = FormPaint
  OnShow = FormShow
  PixelsPerInch = 96
  TextHeight = 13
  object Image1: TImage
    Left = 0
    Top = 0
    Width = 600
    Height = 570
  end
  object StartButton: TButton
    Left = 608
    Top = 568
    Width = 75
    Height = 25
    Caption = 'Start'
    TabOrder = 0
    OnClick = StartButtonClick
  end
  object GraphSelector: TCheckListBox
    Left = 608
    Top = 152
    Width = 185
    Height = 409
    OnClickCheck = GraphSelectorClickCheck
    Flat = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ItemHeight = 13
    ParentFont = False
    TabOrder = 1
  end
  object ProgressBar1: TProgressBar
    Left = 0
    Top = 576
    Width = 601
    Height = 16
    Min = 0
    Max = 1000
    TabOrder = 2
  end
  object EquationSelector: TListBox
    Left = 608
    Top = 0
    Width = 185
    Height = 97
    ItemHeight = 13
    TabOrder = 3
    OnClick = EquationSelectorClick
  end
  object Check_SloMo: TCheckBox
    Left = 608
    Top = 104
    Width = 97
    Height = 17
    Caption = 'Slow Motion'
    TabOrder = 4
  end
  object DumpButton: TButton
    Left = 720
    Top = 568
    Width = 75
    Height = 25
    Caption = 'Show Dump'
    TabOrder = 5
    OnClick = DumpButtonClick
  end
  object Timer1: TTimer
    Interval = 300
    OnTimer = Timer1Timer
    Left = 8
    Top = 536
  end
end
