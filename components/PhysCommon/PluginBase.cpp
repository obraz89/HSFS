///////////////////////////////////////////////////////////////////////////////
// Name:        PluginBase.cpp
// Purpose:     Default behaviour of base classes for plugins
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include <vector>

#include <wx/fileconf.h>

#include <wx/dynlib.h>  //dynamic library support
#include <wx/utils.h>   //wxGetEnv

#include "PluginBase.h"

using namespace hsstab;

//
/////////////////////////////// Constants /////////////////////////////////////
//
static const wxChar strEXE_ENV_VAR[] = _T("HSFLOWStabBin");
static const wxString strDIR_EXE_PLUGINS = _T("");

///////////////////////////////////////////////////////////////////////////////
// class TPluginParam
///////////////////////////////////////////////////////////////////////////////

bool TPluginParam::set_value(int val)
{
	if(type!=ptInt)  return false;
	return _value.Printf(_T("%i"), val) > 0;
}
bool TPluginParam::set_value(double val)
{
	if(type!=ptDouble)  return false;
	return _value.Printf(_T("%G"), val) > 0;
}
bool TPluginParam::set_value(const wxString& aVal)
{
	if(type!=ptString) return false;
	_value = aVal;  return true;
}

//--------

const wxString& TPluginParam::get_raw_value() const
{
	return _value;
}

bool TPluginParam::get_value(int& val) const
{
	if(type!=ptInt)	return false;
	return _value.ToLong((long*)&val);
}
bool TPluginParam::get_value(double& val) const
{
	if(type!=ptDouble) return false;
	return _value.ToDouble(&val);
}
bool TPluginParam::get_value(const wxString*& val) const
{
	if(type!=ptString) return false;

	val = &_value;  return true;
}

//--------

TPluginParam::TPluginParam(){ };

TPluginParam::TPluginParam(int val, const wxString& aDescr/* = wxEmptyString*/)
: description(aDescr), type(ptInt)
{
	set_value(val);
}
TPluginParam::TPluginParam(double val, const wxString& aDescr/* = wxEmptyString*/)
: description(aDescr), type(ptDouble)
{
	set_value(val);
}
TPluginParam::TPluginParam(const wxString& val, const wxString& aDescr/* = wxEmptyString*/)
: description(aDescr), type(ptString)
{
	set_value(val);
}

///////////////////////////////////////////////////////////////////////////////
// class TPluginParamsGroup
///////////////////////////////////////////////////////////////////////////////

TPluginParamsGroup::TPluginParamsGroup(){ };

TPluginParamsGroup::TPluginParamsGroup(const char* aName, 
									   const wxString& aDescr/* = wxEmptyString*/)
: m_description(aDescr)
{
	m_name = wxString::FromAscii(aName);
}
const wxString& TPluginParamsGroup::get_name() const {  return m_name;  }

void TPluginParamsGroup::add(const char* pszParName, double value, 
							 const wxString& aDescr/* = wxEmptyString*/)
{
	wxString parName = wxString::FromAscii(pszParName);
	wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
		_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
		);
	mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
}
void TPluginParamsGroup::add(const char* pszParName, int value, 
							 const wxString& aDescr /* = wxEmptyString*/)
{
	wxString parName = wxString::FromAscii(pszParName);
	wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
		_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
		);
	mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
}
void TPluginParamsGroup::add(const char* pszParName, const wxString& value, 
							 const wxString& aDescr/* = wxEmptyString*/)
{
	wxString parName = wxString::FromAscii(pszParName);
	wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
		_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
		);
	mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
}



///////////////////////////////////////////////////////////////////////////////
// class TPlugin
///////////////////////////////////////////////////////////////////////////////

TPlugin::TPlugin(){  default_settings();  }

void  TPlugin::default_settings()
{
	_mapParamsGrps.clear();
}

TPlugin::~TPlugin(){ ; }


hsstab::TPlugin* hsstab::TPlugin::create_from_dll(const wxString& name)  
{
	static std::vector<wxDllType> loadedPlugs;

	wxString dir;
	if( ! ::wxGetEnv(strEXE_ENV_VAR, &dir) )
		ssuGENTHROW( _("Required environment variable '%s' is not set"), strEXE_ENV_VAR );

	wxFileName plugFN(dir, name, wxDynamicLibrary::GetDllExt()+1/*skip leading dot*/);
	wxLogMessage(_T("loading %s"), plugFN.GetFullPath());
//	plugFN.AppendDir(strDIR_EXE_PLUGINS);

	wxDynamicLibrary lib;
	if( ! lib.Load(plugFN.GetFullPath(), wxDL_VERBATIM|wxDL_NOW|wxDL_GLOBAL) )
		ssuGENTHROW( _("Can't load plugin '%s'"), name.c_str() );
		


	// Check if the same dyn.library has been previously loaded,
	// so here it's not really loaded once more (just ref.counter is increased)
	// Each plugin should have its own address space, they use global variables :(
	wxDllType hndl = lib.GetLibHandle();
	
	for( int i=0; i < loadedPlugs.size(); ++i )
	{
		if( loadedPlugs[i] == hndl )
		{
			ssuGENTHROW( _("Plugin '%s' appears to be loaded twice. Please make a copy of it"), name.c_str() );
		}
	}
	loadedPlugs.push_back( hndl );

	typedef hsstab::TPlugin* (*TpfunPlgIface)();
	TpfunPlgIface get_iface = (TpfunPlgIface)lib.GetSymbol( _T("get_plugin_interface") );
	if( ! get_iface )
		ssuGENTHROW( _("%s: Unable to obtain function 'get_plugin_interface'"), name.c_str() );

	lib.Detach(); // don't unload library

	return get_iface();
}
//-----------------------------------------------------------------------------


hsstab::TPluginParamsGroup& hsstab::TPlugin::get_settings_grp(const char* pszGrpName) 
{
	wxString grpName = wxString::FromAscii(pszGrpName);
	std::map<wxString, TPluginParamsGroup>::iterator it = _mapParamsGrps.find(grpName);
	if( it == _mapParamsGrps.end() )
		ssuGENTHROW( _("%s: Unable to find settings group '%s'"), get_name().c_str(), grpName.c_str() );

	return (*it).second;
}

const hsstab::TPluginParamsGroup& hsstab::TPlugin::get_settings_grp_const(
	const char* pszGrpName) const {

		wxString grpName = wxString::FromAscii(pszGrpName);
		std::map<wxString, TPluginParamsGroup>::const_iterator it = _mapParamsGrps.find(grpName);
		if( it == _mapParamsGrps.end() )
			ssuGENTHROW( _("%s: Unable to find settings group '%s'"), get_name().c_str(), grpName.c_str() );

		return (*it).second;

};
//-----------------------------------------------------------------------------


void hsstab::TPlugin::load_settings(const wxString& fn) 
{
	if( fn.IsEmpty() )  return;

	_paramsFileName = fn;
	wxFileConfig* conf = new wxFileConfig( get_name(), _T("NovA"),
		fn, wxEmptyString, wxCONFIG_USE_RELATIVE_PATH, wxConvUTF8 );

	conf->SetRecordDefaults(); //write defaults to config file

	wxString path = _T("/");
	path += get_name()+_T("/");
	conf->SetPath(path);

	wxString strVal;
	bool isCfgExists = conf->Read(_T("info"), &strVal, get_description());
	if( ! isCfgExists )
	{
		wxLogWarning(
			_("%s: Settings doesn't exist in the file '%s' -- creating defaults..."),
			get_name().c_str(), _paramsFileName.c_str()
		);
	}

	if( ! _spec.IsEmpty() )  path += _spec+_T("/");
	for( std::map<wxString, TPluginParamsGroup>::iterator grp = _mapParamsGrps.begin(); grp != _mapParamsGrps.end(); grp++ )
	{
		conf->SetPath( path + (*grp).first );
		std::map<wxString, TPluginParam>& params = (*grp).second.mapParams;
		for( std::map<wxString, TPluginParam>::iterator prm = params.begin(); prm != params.end(); prm++ )
		{
			bool ok = conf->Read((*prm).first, &strVal, (*prm).second.get_raw_value());
			if(ok) ok = (*prm).second.set_raw_value(strVal);

			if(! ok && isCfgExists)
			{
				wxLogWarning(
					_("%s: Can't read parameter '%s/%s' from the settings file. Using default value %s"),
					get_name().c_str(), (*grp).first.c_str(), (*prm).first.c_str(),
					(*prm).second.get_raw_value().c_str()
				);
			}
		}
	}

	delete conf;
}
//-----------------------------------------------------------------------------


void hsstab::TPlugin::save_settings(const wxString& file) 
{
	const wxString& fn = (file.IsEmpty()) ?_paramsFileName :file;
	if( fn.IsEmpty() )
		ssuGENTHROW( _("%s: Settings file name is not provided"), get_name().c_str() );

	wxFileConfig* conf = new wxFileConfig(get_name(), _T("NovA"), fn, wxEmptyString, wxCONFIG_USE_RELATIVE_PATH);

	wxString path = _T("/");
	path += get_name()+_T("/");

	//Plugin info, just for human reading
	conf->SetPath(path);
	conf->Write(_T("info"), get_description());
	
	if( ! _spec.IsEmpty() )  path += _spec+_T("/");
	for(std::map<wxString, TPluginParamsGroup>::const_iterator grp = _mapParamsGrps.begin(); grp != _mapParamsGrps.end(); grp++)
	{
		conf->SetPath( path + (*grp).first );
		const std::map<wxString, TPluginParam>& params = (*grp).second.mapParams;
		for(std::map<wxString, TPluginParam>::const_iterator prm = params.begin(); prm != params.end(); prm++ )
		{
			if( ! conf->Write( (*prm).first, (*prm).second.get_raw_value() ) )
				wxLogWarning(
					_("%s: Can't write parameter '%s/%s' to the settings file."),
					get_name().c_str(), (*grp).first.c_str(), (*prm).first.c_str()
				);
		}
	}

	delete conf;
}
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class TPluginParamsGroup
///////////////////////////////////////////////////////////////////////////////

const hsstab::TPluginParam&
hsstab::TPluginParamsGroup::get_raw_param(const wxString& parName) const 
{
	std::map<wxString, TPluginParam>::const_iterator it = mapParams.find(parName);
	if( it == mapParams.end() )
		ssuGENTHROW(
			_("Unable to find parameter '%s' in the group '%s'"),
			parName.c_str(), m_name.c_str()
		);

	return (*it).second;
}
//-----------------------------------------------------------------------------

double hsstab::TPluginParamsGroup::get_real_param(const char* pszName) const 
{
 	wxString parName = wxString::FromAscii(pszName);
	const TPluginParam& param = get_raw_param(parName);

	double val;
	if( ! param.get_value(val) )
		ssuGENTHROW(
			_("Parameter '%s' in the group '%s' has value '%s', which is not a real number."),
			parName.c_str(), m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

int hsstab::TPluginParamsGroup::get_int_param(const char* pszName) const 
{
	wxString parName = wxString::FromAscii(pszName);
	const TPluginParam& param = get_raw_param(parName);

	int val;
	if( ! param.get_value(val) )
		ssuGENTHROW(
			_("Parameter '%s' in the group '%s' has value '%s', which is not an integer."),
			parName.c_str(), m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

const wxString& hsstab::TPluginParamsGroup::get_string_param(const char* pszName) const 
{
	wxString parName = wxString::FromAscii(pszName);
	const TPluginParam& par = get_raw_param(parName);

	const wxString* val = NULL;
	if( ! par.get_value(val) )
		ssuGENTHROW(
			_("Parameter '%s' in the group '%s' has value '%s', which is not of string type."),
			parName.c_str(), m_name.c_str(), par.get_raw_value().c_str()
		);

	return *val;
}
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class TPluginParam
///////////////////////////////////////////////////////////////////////////////

bool hsstab::TPluginParam::set_raw_value(const wxString& aVal)
{
	if(type==ptInt)
	{
		long i;
		if( ! aVal.ToLong(&i) )
		{
			wxLogWarning( _("Provided value '%s' is not of required type 'int'"), aVal.c_str() );
			return false;
		}
	}
	else if(type==ptDouble)
	{
		double d;
		if( ! aVal.ToDouble(&d) )
		{
			wxLogWarning( _("Provided value '%s' is not of required type 'double'"), aVal.c_str() );
			return false;
		}
	}

	_value = aVal;  _value.Trim();
	return true;
}
//-----------------------------------------------------------------------------
