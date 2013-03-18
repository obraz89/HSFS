///////////////////////////////////////////////////////////////////////////////
// Name:        PluginBase.h
// Purpose:     Base classes for plugins
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <map>

#include "dll_impexp-phys_common.h"

#include "gen_exception.h"

#pragma warning(disable : 4290) //C++ exception specification ignored
//-----------------------------------------------------------------------------


namespace hsstab
{

enum TPluginType { plgMF=0, plgLS, plgGS, plgWPTrack, plgNUM };

//
// Plugin parameter
// ----------------------------------------------------------------------------
class IMPEXP_PHYSCOMMON TPluginParam
{
private:
	wxString _value;

public:
	enum {ptInt=0, ptDouble, ptString} type;
	wxString description;

public:
	bool set_raw_value(const wxString& aVal);

	bool set_value(int val);
	bool set_value(double val);
	bool set_value(const wxString& aVal);


	const wxString& get_raw_value() const;

	bool get_value(int& val) const;
	bool get_value(double& val) const;
	bool get_value(const wxString*& val) const;


	TPluginParam();
	TPluginParam(int val, const wxString& aDescr = wxEmptyString);
	TPluginParam(double val, const wxString& aDescr = wxEmptyString);
	TPluginParam(const wxString& val, const wxString& aDescr = wxEmptyString);
};
//-----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// class TPluginParamsGroup
//
// Logical group of plugin parameters
// ----------------------------------------------------------------------------
#pragma warning(push)
#pragma warning(disable:1744) //field of class type without a DLL interface used in a class with a DLL interface

class IMPEXP_PHYSCOMMON TPluginParamsGroup
{
	friend class TPlugin;

private:
	wxString m_name;
	wxString m_description;
	std::map<wxString, TPluginParam> mapParams;

	const TPluginParam& get_raw_param(const wxString& parName) const throw(t_GenException);

public:
	TPluginParamsGroup();
	TPluginParamsGroup(const char* aName, const wxString& aDescr = wxEmptyString);
	const wxString& get_name() const;

	void add(const char* pszParName, double value, const wxString& aDescr = wxEmptyString);
	void add(const char* pszParName, int value, const wxString& aDescr = wxEmptyString);
	void add(const char* pszParName, const wxString& value, const wxString& aDescr = wxEmptyString);

	double get_real_param(const char* pszName) const throw(t_GenException);
	int get_int_param(const char* pszName) const throw(t_GenException);
	const wxString& get_string_param(const char* pszName) const throw(t_GenException);
};
#pragma warning(pop)
//-----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// class TPluginCaps
//
// Plugin capabilities (functions)
// ----------------------------------------------------------------------------
class IMPEXP_PHYSCOMMON TPluginCaps
{
	//empty
protected:
	virtual void dummy() = 0;  //to be polymorphic type
};


// ----------------------------------------------------------------------------
// class TPlugin
//
// A base class for all plugin classes.
// ----------------------------------------------------------------------------
#pragma warning(push)
#pragma warning(disable:1744) //field of class type without a DLL interface used in a class with a DLL interface

class IMPEXP_PHYSCOMMON TPlugin
{
protected:
	wxString _spec;  // specialization
	wxString _paramsFileName;
	std::map<wxString, TPluginParamsGroup> _mapParamsGrps;

	virtual void default_settings();

public:
	static TPlugin* create_from_dll(const wxString& name)  throw(t_GenException);

	TPlugin();
	virtual ~TPlugin();

	virtual wxString get_name() const = 0;
	virtual wxString get_description() const = 0;


	virtual void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException)
	{
		_spec = spec;
		load_settings(settingsFN);
	}
	std::map<wxString, TPluginParamsGroup>& get_settings()
	{
		return _mapParamsGrps;
	}
	TPluginParamsGroup& get_settings_grp(const char* pszGrpName) throw(t_GenException);
	const TPluginParamsGroup& get_settings_grp_const(const char* pszGrpName) const throw(t_GenException);

	virtual void load_settings(const wxString& file) throw(t_GenException);
	virtual void save_settings(const wxString& file) throw(t_GenException);

//----- 
	virtual TPluginCaps* get_caps(){ return NULL; }
};

// ----------------------------------------------------------------------------
// class TPlugPhysPart
//
// A base class for implementation of what plugin will do for stability comps.
// ----------------------------------------------------------------------------

class TPlugPhysPart{

public:
	virtual void init(const TPlugin& g_plug)=0;

};

#pragma warning(pop)
//-----------------------------------------------------------------------------

} //namespace hsflow
