///////////////////////////////////////////////////////////////////////////////
// Name:        component.cpp (original : ComponentBase.cpp)
// Purpose:     Default behaviour of base classes for plugins
// Author:      Andrey V. Novikov
// Modified by: Obraz
///////////////////////////////////////////////////////////////////////////////

#include <wx/dynlib.h>  //dynamic library support
#include <wx/utils.h>   //wxGetEnv
#include <wx/fileconf.h>

#include "component.h"

//
/////////////////////////////// Constants /////////////////////////////////////
//


///////////////////////////////////////////////////////////////////////////////
// class t_Component
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------


t_ComponentParamsGroup& t_Component::get_settings_grp(const char* pszGrpName) throw(t_EComponent)
{
	wxString grpName = wxString::FromAscii(pszGrpName);
	std::map<wxString, t_ComponentParamsGroup>::iterator it = _mapParamsGrps.find(grpName);
	if( it == _mapParamsGrps.end() )
//		ssuTHROW( _("%s: Unable to find settings group '%s'"), get_name().c_str(), grpName.c_str() );

	return (*it).second;
}
//-----------------------------------------------------------------------------


void t_Component::load_settings(const wxString& fn) throw(t_EComponent)
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
	for( std::map<wxString, t_ComponentParamsGroup>::iterator grp = _mapParamsGrps.begin(); grp != _mapParamsGrps.end(); grp++ )
	{
		conf->SetPath( path + (*grp).first );
		std::map<wxString, t_ComponentParam>& params = (*grp).second.mapParams;
		for( std::map<wxString, t_ComponentParam>::iterator prm = params.begin(); prm != params.end(); prm++ )
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


void t_Component::save_settings(const wxString& file) throw(t_EComponent)
{
	const wxString& fn = (file.IsEmpty()) ?_paramsFileName :file;
	if( fn.IsEmpty() )
		ssuTHROW( _("%s: Settings file name is not provided"), get_name().c_str() );
		
	wxFileConfig* conf = new wxFileConfig(get_name(), _T("NovA"), fn, wxEmptyString, wxCONFIG_USE_RELATIVE_PATH);

	wxString path = _T("/");
	path += get_name()+_T("/");

	//Component info, just for human reading
	conf->SetPath(path);
	conf->Write(_T("info"), get_description());
	
	if( ! _spec.IsEmpty() )  path += _spec+_T("/");
	for(std::map<wxString, t_ComponentParamsGroup>::const_iterator grp = _mapParamsGrps.begin(); grp != _mapParamsGrps.end(); grp++)
	{
		conf->SetPath( path + (*grp).first );
		const std::map<wxString, t_ComponentParam>& params = (*grp).second.mapParams;
		for(std::map<wxString, t_ComponentParam>::const_iterator prm = params.begin(); prm != params.end(); prm++ )
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
// class t_ComponentParamsGroup
///////////////////////////////////////////////////////////////////////////////

const t_ComponentParam&
t_ComponentParamsGroup::get_raw_param(const wxString& parName) const throw(t_EComponent)
{
	std::map<wxString, t_ComponentParam>::const_iterator it = mapParams.find(parName);
	if( it == mapParams.end() )
		ssuTHROW(
			_("Unable to find parameter '%s' in the group '%s'"),
			parName.c_str(), m_name.c_str()
		);

	return (*it).second;
}
//-----------------------------------------------------------------------------

double t_ComponentParamsGroup::get_real_param(const char* pszName) const throw(t_EComponent)
{
 	wxString parName = wxString::FromAscii(pszName);
	const t_ComponentParam& param = get_raw_param(parName);

	double val;
	if( ! param.get_value(val) )
		ssuTHROW(
			_("Parameter '%s' in the group '%s' has value '%s', which is not a real number."),
			parName.c_str(), m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

int t_ComponentParamsGroup::get_int_param(const char* pszName) const throw(t_EComponent)
{
	wxString parName = wxString::FromAscii(pszName);
	const t_ComponentParam& param = get_raw_param(parName);

	int val;
	if( ! param.get_value(val) )
		ssuTHROW(
			_("Parameter '%s' in the group '%s' has value '%s', which is not an integer."),
			parName.c_str(), m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

/** 
 *  @param pszName   Name of the param to get
 *  @param pRefName  try to dereference parameter and assign pRefName if succeeded
 */
const wxString& t_ComponentParamsGroup::get_string_param(const char* pszName, wxString* pRefName) const throw(t_EComponent)
{
	wxString parName = wxString::FromAscii(pszName);
	const t_ComponentParam& par = get_raw_param(parName);

	const wxString* val = NULL;
	if( ! par.get_value(val) )   goto err;

	if(pRefName)  // try to dereference
	{
		if( val->StartsWith(_T("@"), pRefName) )  // dereference
		{
			try
			{
				const t_ComponentParam& par2 = get_raw_param(*pRefName);
				if( ! par2.get_value(val) )  goto err;
			}
			catch(const t_EComponent& e)
			{
				*pRefName = wxEmptyString;
			}
		}
	}

	return *val;

err:
	ssuTHROW(
		_("Parameter '%s' in the group '%s' has value '%s', which is not of string type."),
		parName.c_str(), m_name.c_str(), par.get_raw_value().c_str()
	);
}
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class t_ComponentParam
///////////////////////////////////////////////////////////////////////////////

bool t_ComponentParam::set_raw_value(const wxString& aVal)
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
