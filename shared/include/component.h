///////////////////////////////////////////////////////////////////////////////
// Name:        component.h (original : PluginBase.)
// Purpose:     Base classes for plugins/components
// Author:      Andrey V. Novikov
// Modified by: Obraz
// Note: syntax modified
///////////////////////////////////////////////////////////////////////////////
#ifndef __MY_COMPONENT
#define __MY_COMPONENT
#include "stdafx.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "wx/fileconf.h"
#include "gen_exception.h"

class t_BaseParam{
public:
	enum t_Type{Int, Dbl, Str};
protected:
	wxString _name;
public:
	wxString description;
	t_BaseParam(wxString name, wxString desc=wxEmptyString):_name(name), description(desc){};
	const wxString& name() const{return _name;};
	virtual t_Type type()=0;
};

// Component param is an interface (or better say proxy) to
// a certain really existing param.
// this provides link to Component param hierarchy
/*
template<typename T> class t_ComponentParam : public t_BaseParam{
private:
	T* _pVal;
	T _defaultVal;
public:
	inline void set_value(T val){*_pVal=val;};
	inline T get_value(){return *_pVal};
	inline void set_default(T val){_defaultVal = val;};
	inline T get_default(){return _defaultVal;};
	t_ComponentParam(t_BaseParam::t_Type type, T& val, wxString name,wxString desc=wxEmptyString):
		t_BaseParam(type, name), _pVal(&val), description(desc){};
		
	void set_value_str(wxString strVal){
			std::string stdstr(strVal.c_str());
			T val;
			try{
				val =  boost::lexical_cast<T, std::string>(stdstr);
			}
			catch(boost::bad_lexical_cast){
				wxLogError(_("Unsupported conversion in ComponentParam initialization"));
				//std::cerr<<"Error parsing input string:\n\t"<<s_out.str()<<std::endl;
				// abort ...
			}
			*_pVal = val;
	};
	
	//exceptions
};
*/
class t_Enum{
	friend class t_CompParamInt;
protected:
	std::map<int, wxString> _mapVals;
	virtual void _init_map_vals()=0;
	int _curVal;
	int* _get_val_addr(){return &_curVal;};
public:
	virtual void set_value(int val){
		// enum)))
		if(_mapVals.find(val)==_mapVals.end()) return;
		_curVal=val;
	};
	virtual int get_value(){return _curVal;};
	virtual bool operator==(int val) const{return _curVal==val;};
	virtual void operator=(const int& val){set_value(val);};
};

class t_CompParamInt : public t_BaseParam{
private:
	int* _pVal;
	int _defaultVal;
public:
	inline void set_value(int val){*_pVal=val;};
	inline int get_value(){return *_pVal;};
	inline void set_default(int val){_defaultVal = val;};
	inline int get_default(){return _defaultVal;};
	t_Type type(){return t_Type::Int;};
	t_CompParamInt(int& val, wxString name,wxString desc=wxEmptyString):
		t_BaseParam(name, desc), _pVal(&val){};
	// reference to enum
	t_CompParamInt(t_Enum& enum_val, wxString name,wxString desc=wxEmptyString):
			t_BaseParam(name, desc), _pVal(enum_val._get_val_addr()){};
};

class t_CompParamDbl : public t_BaseParam{
private:
	double* _pVal;
	double _defaultVal;
public:
	inline void set_value(double val){*_pVal=val;};
	inline double get_value(){return *_pVal;};
	inline void set_default(double val){_defaultVal = val;};
	inline double get_default(){return _defaultVal;};
	t_Type type(){return t_Type::Dbl;};
	t_CompParamDbl(double& val, wxString name,wxString desc=wxEmptyString):
	t_BaseParam(name, desc), _pVal(&val){};
};

class t_CompParamStr : public t_BaseParam{
private:
	wxString* _pVal;
	wxString _defaultVal;
public:
	inline void set_value(wxString val){*_pVal=val;};
	inline wxString get_value(){return *_pVal;};
	inline void set_default(wxString val){_defaultVal = val;};
	inline wxString get_default(){return _defaultVal;};
	t_Type type(){return t_Type::Str;};
	t_CompParamStr(wxString& val, wxString name,wxString desc=wxEmptyString):
	t_BaseParam(name, desc), _pVal(&val){};
};

/*
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
*/

class t_ComponentParamsGroup{
protected:
	const wxString ConfigDomain;
	std::map<wxString, t_BaseParam*> _mapParams;
	virtual void _init_params_map()=0;
	void _add_param(t_BaseParam* pParam);
	wxFileConfig _get_config_handle(wxString configfile);
	virtual void _load_via_params(wxFileConfig& handle);
	void _clear_map();
public:
	t_ComponentParamsGroup(wxString config_domain):ConfigDomain(config_domain){};
	virtual ~t_ComponentParamsGroup(){
		_clear_map();
	};
	virtual void load_direct(wxString configfile)=0;
	virtual void load_via_params(wxString configfile)=0;
	virtual void save(wxString configfile)=0;
};
/*
class t_ComponentParamsGroup
{
	friend class t_Component;

private:
	wxString m_name;
	wxString m_description;
	std::map<wxString, t_ComponentParam> mapParams;

	const t_ComponentParam& get_raw_param(const wxString& parName) const throw(t_GenException);

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

	double get_real_param(const char* pszName) const throw(t_GenException);
	int get_int_param(const char* pszName) const throw(t_GenException);
	const wxString& get_string_param(const char* pszName, wxString* pRefName = NULL) const throw(t_GenException);
};

*/
class t_Component{
protected:
	std::map<wxString, t_ComponentParamsGroup*> _mapParamsGrps;
	wxString _paramsFileName;
	wxString _name;
	wxString _spec;
	virtual void _init_params_grps()=0;
	void _add_params_group(wxString name, t_ComponentParamsGroup&);
public:
	virtual ~t_Component(){ ; }
	const wxString& name() const{return _name;};
	const wxString& spec() const{return _spec;};


	t_Component(wxString settingsFN, wxString name,wxString spec=wxEmptyString) throw(t_GenException)
	{
		_paramsFileName=settingsFN;
		_name = name;
		_spec = spec;
		//	load_settings(settingsFN);
	};
	/*
	std::map<wxString, t_ComponentParamsGroup>& get_settings()
	{
		return _mapParamsGrps;
	}*/
	/*
	t_ComponentParamsGroup& get_settings_grp(const char* pszGrpName) throw(t_GenException);

	virtual wxString get_name() const = 0;
	virtual wxString get_description() const = 0;
	*/
	virtual void initialize(const wxString& file)=0;
	virtual void load_settings(const wxString& file) throw(t_GenException);
	virtual void save_settings(const wxString& file) throw(t_GenException);
};

#endif // __MY_COMPONENT