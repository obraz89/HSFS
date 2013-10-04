///////////////////////////////////////////////////////////////////////////////
// Name:        io.h
// Purpose:     Input/Output functions
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

bool initField();
bool loadField(const wxString& fileName, bool isPrev);

bool saveField(const wxString& fileName, const wxString& gridFileName);
void saveField_legacy(const char* aFileName);
//-----------------------------------------------------------------------------
