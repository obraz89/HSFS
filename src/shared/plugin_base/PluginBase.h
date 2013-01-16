///////////////////////////////////////////////////////////////////////////////
// Name:        PluginBase.h
// Purpose:     Base classes for plugins
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
///////////////////////////////////////////////////////////////////////////////


// Plugin is a configurable unit of the program

#pragma once

#include <map>

#include "dll_import-export.h"

//#pragma warning(disable : 4290) //C++ exception specification ignored
//-----------------------------------------------------------------------------


#define genTHROW(...)  \
	throw EPlugin( wxString::Format(__VA_ARGS__), __TFILE__, __LINE__ )
//-----------------------------------------------------------------------------




	//-----------------------------------------------------------------------------
	// class EPlugin
	//
	// A base class for all plugin exceptions
	//-----------------------------------------------------------------------------
	class DLLIMPEXP EPlugin
	{
	protected:
		wxString    _what;
		wxString    _file;
		int         _line;

	public:
		EPlugin(const wxString& what, const wxChar* szFile,  const int line)
			:_what(what), _file(szFile), _line(line){};

		wxString what() const  {  return _what;  }
		wxString file() const  {  return _file;  }
		int line() const       {  return _line;  }
		wxString what_detailed() const
		{
			return _what + _(". In file: ")+_file + _(", line: ")+wxString::Format(_T("%d"), _line);
		}
	};
	//-----------------------------------------------------------------------------


	// ----------------------------------------------------------------------------
	// class TPluginParam
	//
	// Plugin parameter
	// ----------------------------------------------------------------------------
	class DLLIMPEXP TPluginParam
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

		const wxString& get_raw_value() const;

		bool get_value(int& val) const;
		bool get_value(double& val) const;
		bool get_value(const wxString*& val) const;
		//--------

		TPluginParam(){ };
		TPluginParam(int val, const wxString& aDescr = wxEmptyString)
			: description(aDescr), type(ptInt)
		{
			set_value(val);
		}
		TPluginParam(double val, const wxString& aDescr = wxEmptyString)
			: description(aDescr), type(ptDouble)
		{
			set_value(val);
		}
		TPluginParam(const wxString& val, const wxString& aDescr = wxEmptyString)
			: description(aDescr), type(ptString)
		{
			set_value(val);
		}
	};
	//-----------------------------------------------------------------------------


	// ----------------------------------------------------------------------------
	// class TPluginParamsGroup
	//
	// Logical group of plugin parameters
	// ----------------------------------------------------------------------------
#pragma warning(push)
#pragma warning(disable:1744) //field of class type without a DLL interface used in a class with a DLL interface

	class DLLIMPEXP TPluginParamsGroup
	{
		friend class TPlugin;

	private:
		wxString m_name;
		wxString m_description;
		std::map<wxString, TPluginParam> mapParams;

		const TPluginParam& get_raw_param(const wxString& parName) const throw(EPlugin);

	public:
		TPluginParamsGroup(){ };
		TPluginParamsGroup(const char* aName, const wxString& aDescr = wxEmptyString)
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
			mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
		}
		void add(const char* pszParName, int value, const wxString& aDescr = wxEmptyString)
		{
			wxString parName = wxString::FromAscii(pszParName);
			wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
				_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
				);
			mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
		}
		void add(const char* pszParName, const wxString& value, const wxString& aDescr = wxEmptyString)
		{
			wxString parName = wxString::FromAscii(pszParName);
			wxASSERT_MSG( mapParams.find(parName)==mapParams.end(),
				_("Parameter '")+parName+_("' already exists in the group '")+m_name+_("'.")
				);
			mapParams.insert( std::make_pair(parName, TPluginParam(value, aDescr)) );
		}

		double get_real_param(const char* pszName) const throw(EPlugin);
		int get_int_param(const char* pszName) const throw(EPlugin);
		const wxString& get_string_param(const char* pszName) const throw(EPlugin);
	};
#pragma warning(pop)
	//-----------------------------------------------------------------------------


	// ----------------------------------------------------------------------------
	// class TPluginCaps
	//
	// Plugin capabilities (functions)
	// ----------------------------------------------------------------------------
	class DLLIMPEXP TPluginCaps
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

	class DLLIMPEXP TPlugin
	{
	protected:
		wxString _spec;  // specialization
		wxString _paramsFileName;
		std::map<wxString, TPluginParamsGroup> _mapParamsGrps;

		virtual void default_settings()
		{
			_mapParamsGrps.clear();
		}

	public:
		static TPlugin* create_from_dll(const wxString& name)  throw(EPlugin);

		TPlugin(){  default_settings();  }
		virtual ~TPlugin(){ ; }

		virtual wxString get_name() const = 0;
		virtual wxString get_description() const = 0;


		virtual void init(const wxString& settingsFN, const wxString& spec) throw(EPlugin)
		{
			_spec = spec;
			load_settings(settingsFN);
		}
		std::map<wxString, TPluginParamsGroup>& get_settings()
		{
			return _mapParamsGrps;
		}
		TPluginParamsGroup& get_settings_grp(const char* pszGrpName) throw(EPlugin);
		virtual void load_settings(const wxString& file) throw(EPlugin);
		virtual void save_settings(const wxString& file) throw(EPlugin);

		//----- 
		virtual TPluginCaps* get_caps(){ return NULL; }
	};

#pragma warning(pop)
