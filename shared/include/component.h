///////////////////////////////////////////////////////////////////////////////
// Name:        component.h (original : PluginBase.)
// Purpose:     Base classes for plugins/components
// Author:      Andrey V. Novikov
// Modified by: Obraz
// Note: synthax modified
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#define ssuTHROW(...)  \
	throw t_EComponent( wxString::Format(__VA_ARGS__), __TFILE__, __LINE__ )
//------

class t_EComponent
{
protected:
	wxString    _what;
	wxString    _file;
	int         _line;

public:
	t_EComponent(const wxString& what, const wxChar* szFile,  const int line)
		: _what(what), _file(szFile), _line(line) {   }

	wxString what() const  {  return _what;  }
	wxString file() const  {  return _file;  }
	int line() const       {  return _line;  }
	wxString what_detailed() const
	{
		return _what + _(". In file: ")+_file + _(", line: ")+wxString::Format(_T("%d"), _line);
	}
};

class t_ComponentParam
{
private:
	wxString _value;

public:
	enum {ptInt=0, ptDouble, ptString} type;
	wxString description;

public:
	bool set_raw_value(const wxString& aVal);

	bool set_value(int val)
	{
		if(type!=ptInt)  return false;
		return _value.Printf(_T("%i"), val) > 0;
	}
	bool set_value(double val)
	{
		if(type!=ptDouble)  return false;
		return _value.Printf(_T("%G"), val) > 0;
	}
	bool set_value(const wxString& aVal)
	{
		if(type!=ptString) return false;
		_value = aVal;  return true;
	}

	//--------

	const wxString& get_raw_value() const
	{
		return _value;
	}

	bool get_value(int& val) const
	{
		if(type!=ptInt)	return false;
		return _value.ToLong((long*)&val);
	}
	bool get_value(double& val) const
	{
		if(type!=ptDouble) return false;
		return _value.ToDouble(&val);
	}
	bool get_value(const wxString*& val) const
	{
		if(type!=ptString) return false;

		val = &_value;  return true;
	}

	//--------

	t_ComponentParam(){ };
	t_ComponentParam(int val, const wxString& aDescr = wxEmptyString)
		: description(aDescr), type(ptInt)
	{
		set_value(val);
	}
	t_ComponentParam(double val, const wxString& aDescr = wxEmptyString)
		: description(aDescr), type(ptDouble)
	{
		set_value(val);
	}
	t_ComponentParam(const wxString& val, const wxString& aDescr = wxEmptyString)
		: description(aDescr), type(ptString)
	{
		set_value(val);
	}
};
//-----------------------------------------------------------------------------

class t_ComponentParamsGroup
{
	friend class t_Component;

private:
	wxString m_name;
	wxString m_description;
	std::map<wxString, t_ComponentParam> mapParams;

	const t_ComponentParam& get_raw_param(const wxString& parName) const throw(t_EComponent);

public:
	t_ComponentParamsGroup(){ };
	t_ComponentParamsGroup(const char* aName, const wxString& aDescr = wxEmptyString)
		: m_description(aDescr)
	{
		m_name = wxString::FromAscii(aName);
	}
	const wxString& get_name() const {  return m_name;  }

	void add(const char* pszParName, double value, const wxString& aDescr = wxEmptyString)
	{
		wxString parName = wxString::FromAscii(pszParName);
		wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
			_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
			);
		mapParams.insert( std::make_pair(parName, t_ComponentParam(value, aDescr)) );
	}
	void add(const char* pszParName, int value, const wxString& aDescr = wxEmptyString)
	{
		wxString parName = wxString::FromAscii(pszParName);
		wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
			_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
			);
		mapParams.insert( std::make_pair(parName, t_ComponentParam(value, aDescr)) );
	}
	void add(const char* pszParName, const wxString& value, const wxString& aDescr = wxEmptyString)
	{
		wxString parName = wxString::FromAscii(pszParName);
		wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
			_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
			);
		mapParams.insert( std::make_pair(parName, t_ComponentParam(value, aDescr)) );
	}

	double get_real_param(const char* pszName) const throw(t_EComponent);
	int get_int_param(const char* pszName) const throw(t_EComponent);
	const wxString& get_string_param(const char* pszName, wxString* pRefName = NULL) const throw(t_EComponent);
};


class t_Component{
	std::map<wxString, t_ComponentParamsGroup> _mapParamsGrps;
	wxString _spec;  // specialization
	wxString _paramsFileName;

	virtual void default_settings()
	{
		_mapParamsGrps.clear();
	}

	t_Component(){  default_settings();  }
	virtual ~t_Component(){ ; }


	virtual void init(const wxString& settingsFN, const wxString& spec) throw(t_EComponent)
	{
		_spec = spec;
		load_settings(settingsFN);
	}
	std::map<wxString, t_ComponentParamsGroup>& get_settings()
	{
		return _mapParamsGrps;
	}
	t_ComponentParamsGroup& get_settings_grp(const char* pszGrpName) throw(t_EComponent);
	virtual void load_settings(const wxString& file) throw(t_EComponent);
	virtual void save_settings(const wxString& file) throw(t_EComponent);

	virtual wxString get_name() const = 0;
	virtual wxString get_description() const = 0;

};