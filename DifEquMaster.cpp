//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
//---------------------------------------------------------------------------
USEFORM("MainForm.cpp", Form1);
USEFORM("Equation01.cpp", DummyForm);
USEFORM("Equation02.cpp", DummyForm2);
USEFORM("OpticParaBasicEqu.cpp", FormOpticBasic);
USEFORM("Equ25DEl01.cpp", FormEqu25DEl01);
USEFORM("EquBaseDump.cpp", DumpForm);
USEFORM("Dump25DElForm.cpp", Dump25DEl);
USEFORM("EquOptic2DMasterForm.cpp", EquOptic2DForm);
USEFORM("Dump2DOpticsForm.cpp", Dump2DOpticsForm);
USEFORM("PSOFormUnit.cpp", PSOForm);
USEFORM("PSOSolveMasterForm.cpp", PSOSolverForm);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
        try
        {
                 Application->Initialize();
                 Application->CreateForm(__classid(TForm1), &Form1);
                 Application->CreateForm(__classid(TDummyForm), &DummyForm);
                 Application->CreateForm(__classid(TDummyForm2), &DummyForm2);
                 Application->CreateForm(__classid(TFormOpticBasic), &FormOpticBasic);
                 Application->CreateForm(__classid(TFormEqu25DEl01), &FormEqu25DEl01);
                 Application->CreateForm(__classid(TDumpForm), &DumpForm);
                 Application->CreateForm(__classid(TDump25DEl), &Dump25DEl);
                 Application->CreateForm(__classid(TEquOptic2DForm), &EquOptic2DForm);
                 Application->CreateForm(__classid(TPSOForm), &PSOForm);
                 Application->CreateForm(__classid(TPSOSolverForm), &PSOSolverForm);
                 Application->Run();
        }
        catch (Exception &exception)
        {
                 Application->ShowException(&exception);
        }
        catch (...)
        {
                 try
                 {
                         throw Exception("");
                 }
                 catch (Exception &exception)
                 {
                         Application->ShowException(&exception);
                 }
        }
        return 0;
}
//---------------------------------------------------------------------------
