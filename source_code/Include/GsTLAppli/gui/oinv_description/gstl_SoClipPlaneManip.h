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

#ifndef __GSTLAPPLI_OINV_GSTLSOCLIPPLANEMANIP_H__ 
#define __GSTLAPPLI_OINV_GSTLSOCLIPPLANEMANIP_H__ 


#include <GsTLAppli/gui/common.h>


#ifdef WIN32
#include <SoWinLeaveScope.h>
#endif


#include <Inventor/manips/SoClipPlaneManip.h> 
#include <Inventor/fields/SoSFBool.h> 
#include <Inventor/fields/SoSFInt32.h> 
 
class SoAction; 
class SoGLRenderAction; 
 

class GUI_DECL GsTL_SoClipPlaneManip : public SoClipPlaneManip { 
 
  SO_NODE_HEADER(GsTL_SoClipPlaneManip); 
   
 public: 
 
  // Initializes the class 
  static void initClass(); 
 
  // default constructor 
  GsTL_SoClipPlaneManip(); 
  virtual ~GsTL_SoClipPlaneManip() {} 
 
}; 


#ifdef WIN32
#include <SoWinEnterScope.h>
#endif

 
#endif 
