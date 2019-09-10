/**********************************************************************
** Author: Nicolas Remy
** Copyright (C) 2002-2004 The Board of Trustees of the Leland Stanford Junior
**   University
** All rights reserved.
**
** This file is part of the "gui" module of the Geostatistical Earth
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

#ifndef __GSTLAPPLI_GUI_FIRST_SCR_H__
#define __GSTLAPPLI_GUI_FIRST_SCR_H__

#include <GsTLAppli/gui/ui_first_scr_base.h>

#include <GsTLAppli/extra/qtplugins/selectors.h>
#include <GsTLAppli/appli/project.h>
#include <GsTLAppli/math/gstlvector.h>
#include <GsTLAppli/math/discrete_function.h>
#include <GsTLAppli/utils/string_manipulation.h>
#include <qstring.h>
#include <qdom.h>
#include <qgroupbox.h>

#include "common.h"	

#include <vector>	

using namespace String_Op;

class FileChooser;

class /*GUI_DECL*/ First_scr : public QWidget, public Ui::First_scr_base {
  
  Q_OBJECT

public:
  First_scr(){}
  First_scr( GsTL_project* project, 
             QWidget* parent = 0, const char* name = 0 );
  std::string grid_name() const ;
  std::string head_property() const;
  std::string tail_property() const;

  bool skip_variogram_computation();
  std::vector<Discrete_function>& experimental_variograms() {
    return experimental_variograms_;
  }
  QStringList& parameters() { return parameters_list_; }
  std::vector< GsTLVector<double> >& directions() { return directions_; }
  std::vector< std::vector<int> >& pairs() { return pairs_; }

public:// 
  bool load_experimental_variograms(char * filename0 = 0);


signals:
  void parameterChanged(); 


private:
  QString extract_title( QDomElement& elem, int id );
  bool extract_direction( GsTLVector<double>& dir,
                          QDomElement& elem, int id );
  bool extract_values( std::vector<double>& values, const QString& tagname,
                       QDomElement& elem, int id );

private:
  GridSelectorBasic *grid_;
  SinglePropertySelector *head_prop_;
  SinglePropertySelector *tail_prop_;
  FileChooser* filechooser_;

  std::vector<Discrete_function> experimental_variograms_;
  QStringList parameters_list_;
  std::vector< GsTLVector<double> > directions_;
  std::vector< std::vector<int> > pairs_;
};






#endif
