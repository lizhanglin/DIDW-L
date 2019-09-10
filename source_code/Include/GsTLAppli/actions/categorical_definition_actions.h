/**********************************************************************
** Author: Alexandre Boucher
** Copyright (C) 2002-2004 The Board of Trustees of the Leland Stanford Junior
**   University
** All rights reserved.
**
** This file is part of the "actions" module of the Geostatistical Earth
** Modeling Software (GEMS)
**
** This file may be distributed and/or modified under the terms of the 
** license defined by the Stanford Center for Reservoir Forecasting and 
** appearing in the file LICENSE.XFREE included in the packaging of this file.
**
** This file may be distributed and/or modified under the terms of the
** GNU General Public License version 2 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
** See http://www.gnu.org/copyleft/gpl.html for GPL licensing information.
**
** Contact the Stanford Center for Reservoir Forecasting, Stanford University
** if any conditions of this licensing are not clear to you.
**
**********************************************************************/

#ifndef __GSTLAPPLI_CATEGORICAL_DEFINITION_ACTIONS_H__ 
#define __GSTLAPPLI_CATEGORICAL_DEFINITION_ACTIONS_H__ 

#include <GsTLAppli/actions/common.h>
#include <GsTLAppli/actions/action.h> 


#include <qstring.h>
 
#include <string> 
#include <set>

class GsTL_project; 
class Error_messages_handler;

class ACTIONS_DECL New_categorical_definition : public Action { 
 public: 
  static Named_interface* create_new_interface( std::string& ); 
 
 public: 
  New_categorical_definition() {} 
  virtual ~New_categorical_definition() {} 
   
  /** The parameters for this function are: 
   * - the name of the new grid 
   * - its dimensions nx, ny, nz 
   */ 
  virtual bool init( std::string& parameters, GsTL_project* proj,
                     Error_messages_handler* errors ); 
  virtual bool exec(); 
 
 private: 

}; 


class ACTIONS_DECL Assign_categorical_definition : public Action { 
 public: 
  static Named_interface* create_new_interface( std::string& ); 
 
 public: 
  Assign_categorical_definition() {} 
  virtual ~Assign_categorical_definition() {} 
   
  /** The parameters for this function are: 
   * - the name of the new grid 
   * - its dimensions nx, ny, nz 
   */ 
  virtual bool init( std::string& parameters, GsTL_project* proj,
                     Error_messages_handler* errors ); 
  virtual bool exec(); 
 
 private: 

  GsTL_project* proj_; 
  Error_messages_handler* errors_;
}; 



#endif
