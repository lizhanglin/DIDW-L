/********************************************************************************
** Form generated from reading UI file 'point_set_base.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_POINT_SET_BASE_H
#define UI_POINT_SET_BASE_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Point_set_base
{
public:
    QVBoxLayout *vboxLayout;
    QHBoxLayout *hboxLayout;
    QSpacerItem *spacer12;
    QPushButton *load_button_;
    QPushButton *save_button_;
    QSpacerItem *spacer11;
    QSpacerItem *spacer10_3_2;
    QHBoxLayout *hboxLayout1;
    QLabel *textLabel1;
    QFrame *line1;
    QSpacerItem *spacer10_3;
    QHBoxLayout *hboxLayout2;
    QSpacerItem *spacer7;
    QGridLayout *gridLayout;
    QLineEdit *lag_tol_;
    QLabel *TextLabel4;
    QSpinBox *num_lags_;
    QLabel *TextLabel1_2;
    QLabel *TextLabel1;
    QLineEdit *lag_sep_;
    QSpacerItem *Spacer9;
    QLabel *pixmapLabel6;
    QSpacerItem *spacer18;
    QSpacerItem *Spacer10;
    QHBoxLayout *hboxLayout3;
    QLabel *textLabel1_2;
    QFrame *line1_2;
    QSpacerItem *spacer10_2;
    QHBoxLayout *hboxLayout4;
    QSpacerItem *spacer14;
    QLabel *textLabel1_4;
    QSpinBox *directions_count_;
    QSpacerItem *spacer15;
    QLabel *pixmapLabel5;
    QSpacerItem *spacer17;
    QLabel *textLabel1_5;
    QSpacerItem *spacer16;
    QGroupBox *table_frame_;
    QSpacerItem *Spacer10_2;
    QGroupBox *standardize_options_box_;
    QVBoxLayout *vboxLayout1;
    QHBoxLayout *hboxLayout5;
    QLabel *textLabel1_3;
    QFrame *line1_3;
    QSpacerItem *spacer10;
    QHBoxLayout *hboxLayout6;
    QSpacerItem *spacer8_2;
    QCheckBox *standardize_checkbox_;
    QSpacerItem *spacer8;
    QSpacerItem *spacer5;

    void setupUi(QWidget *Point_set_base)
    {
        if (Point_set_base->objectName().isEmpty())
            Point_set_base->setObjectName(QString::fromUtf8("Point_set_base"));
        Point_set_base->resize(637, 656);
        vboxLayout = new QVBoxLayout(Point_set_base);
        vboxLayout->setSpacing(6);
        vboxLayout->setContentsMargins(14, 14, 14, 14);
        vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
        hboxLayout = new QHBoxLayout();
        hboxLayout->setSpacing(6);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        spacer12 = new QSpacerItem(200, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacer12);

        load_button_ = new QPushButton(Point_set_base);
        load_button_->setObjectName(QString::fromUtf8("load_button_"));

        hboxLayout->addWidget(load_button_);

        save_button_ = new QPushButton(Point_set_base);
        save_button_->setObjectName(QString::fromUtf8("save_button_"));

        hboxLayout->addWidget(save_button_);

        spacer11 = new QSpacerItem(31, 20, QSizePolicy::Minimum, QSizePolicy::Minimum);

        hboxLayout->addItem(spacer11);


        vboxLayout->addLayout(hboxLayout);

        spacer10_3_2 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(spacer10_3_2);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setSpacing(6);
        hboxLayout1->setContentsMargins(0, 0, 0, 0);
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        textLabel1 = new QLabel(Point_set_base);
        textLabel1->setObjectName(QString::fromUtf8("textLabel1"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(textLabel1->sizePolicy().hasHeightForWidth());
        textLabel1->setSizePolicy(sizePolicy);
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        textLabel1->setFont(font);
        textLabel1->setWordWrap(false);

        hboxLayout1->addWidget(textLabel1);

        line1 = new QFrame(Point_set_base);
        line1->setObjectName(QString::fromUtf8("line1"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(line1->sizePolicy().hasHeightForWidth());
        line1->setSizePolicy(sizePolicy1);
        line1->setFrameShape(QFrame::HLine);
        line1->setFrameShadow(QFrame::Raised);

        hboxLayout1->addWidget(line1);


        vboxLayout->addLayout(hboxLayout1);

        spacer10_3 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(spacer10_3);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setSpacing(6);
        hboxLayout2->setContentsMargins(0, 0, 0, 0);
        hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
        spacer7 = new QSpacerItem(120, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        hboxLayout2->addItem(spacer7);

        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(0, 0, 0, 0);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        lag_tol_ = new QLineEdit(Point_set_base);
        lag_tol_->setObjectName(QString::fromUtf8("lag_tol_"));
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(lag_tol_->sizePolicy().hasHeightForWidth());
        lag_tol_->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(lag_tol_, 2, 1, 1, 1);

        TextLabel4 = new QLabel(Point_set_base);
        TextLabel4->setObjectName(QString::fromUtf8("TextLabel4"));
        TextLabel4->setWordWrap(false);

        gridLayout->addWidget(TextLabel4, 2, 0, 1, 1);

        num_lags_ = new QSpinBox(Point_set_base);
        num_lags_->setObjectName(QString::fromUtf8("num_lags_"));
        num_lags_->setMinimum(1);
        num_lags_->setMaximum(500);
        num_lags_->setSingleStep(1);

        gridLayout->addWidget(num_lags_, 0, 1, 1, 1);

        TextLabel1_2 = new QLabel(Point_set_base);
        TextLabel1_2->setObjectName(QString::fromUtf8("TextLabel1_2"));
        TextLabel1_2->setWordWrap(false);

        gridLayout->addWidget(TextLabel1_2, 1, 0, 1, 1);

        TextLabel1 = new QLabel(Point_set_base);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setWordWrap(false);

        gridLayout->addWidget(TextLabel1, 0, 0, 1, 1);

        lag_sep_ = new QLineEdit(Point_set_base);
        lag_sep_->setObjectName(QString::fromUtf8("lag_sep_"));
        sizePolicy2.setHeightForWidth(lag_sep_->sizePolicy().hasHeightForWidth());
        lag_sep_->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(lag_sep_, 1, 1, 1, 1);


        hboxLayout2->addLayout(gridLayout);

        Spacer9 = new QSpacerItem(160, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout2->addItem(Spacer9);

        pixmapLabel6 = new QLabel(Point_set_base);
        pixmapLabel6->setObjectName(QString::fromUtf8("pixmapLabel6"));
        pixmapLabel6->setPixmap(QPixmap(QString::fromUtf8("image0")));
        pixmapLabel6->setScaledContents(true);
        pixmapLabel6->setWordWrap(false);

        hboxLayout2->addWidget(pixmapLabel6);

        spacer18 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        hboxLayout2->addItem(spacer18);


        vboxLayout->addLayout(hboxLayout2);

        Spacer10 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(Spacer10);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setSpacing(6);
        hboxLayout3->setContentsMargins(0, 0, 0, 0);
        hboxLayout3->setObjectName(QString::fromUtf8("hboxLayout3"));
        textLabel1_2 = new QLabel(Point_set_base);
        textLabel1_2->setObjectName(QString::fromUtf8("textLabel1_2"));
        sizePolicy.setHeightForWidth(textLabel1_2->sizePolicy().hasHeightForWidth());
        textLabel1_2->setSizePolicy(sizePolicy);
        textLabel1_2->setFont(font);
        textLabel1_2->setWordWrap(false);

        hboxLayout3->addWidget(textLabel1_2);

        line1_2 = new QFrame(Point_set_base);
        line1_2->setObjectName(QString::fromUtf8("line1_2"));
        sizePolicy1.setHeightForWidth(line1_2->sizePolicy().hasHeightForWidth());
        line1_2->setSizePolicy(sizePolicy1);
        line1_2->setFrameShape(QFrame::HLine);
        line1_2->setFrameShadow(QFrame::Raised);

        hboxLayout3->addWidget(line1_2);


        vboxLayout->addLayout(hboxLayout3);

        spacer10_2 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(spacer10_2);

        hboxLayout4 = new QHBoxLayout();
        hboxLayout4->setSpacing(6);
        hboxLayout4->setContentsMargins(0, 0, 0, 0);
        hboxLayout4->setObjectName(QString::fromUtf8("hboxLayout4"));
        spacer14 = new QSpacerItem(120, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        hboxLayout4->addItem(spacer14);

        textLabel1_4 = new QLabel(Point_set_base);
        textLabel1_4->setObjectName(QString::fromUtf8("textLabel1_4"));
        textLabel1_4->setWordWrap(false);

        hboxLayout4->addWidget(textLabel1_4);

        directions_count_ = new QSpinBox(Point_set_base);
        directions_count_->setObjectName(QString::fromUtf8("directions_count_"));
        directions_count_->setButtonSymbols(QAbstractSpinBox::UpDownArrows);
        directions_count_->setMinimum(1);
        directions_count_->setMaximum(50);

        hboxLayout4->addWidget(directions_count_);

        spacer15 = new QSpacerItem(78, 21, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout4->addItem(spacer15);

        pixmapLabel5 = new QLabel(Point_set_base);
        pixmapLabel5->setObjectName(QString::fromUtf8("pixmapLabel5"));
        pixmapLabel5->setPixmap(QPixmap(QString::fromUtf8("image1")));
        pixmapLabel5->setScaledContents(true);
        pixmapLabel5->setWordWrap(false);

        hboxLayout4->addWidget(pixmapLabel5);

        spacer17 = new QSpacerItem(16, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        hboxLayout4->addItem(spacer17);

        textLabel1_5 = new QLabel(Point_set_base);
        textLabel1_5->setObjectName(QString::fromUtf8("textLabel1_5"));
        textLabel1_5->setMinimumSize(QSize(150, 0));
        textLabel1_5->setWordWrap(false);

        hboxLayout4->addWidget(textLabel1_5);


        vboxLayout->addLayout(hboxLayout4);

        spacer16 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(spacer16);

        table_frame_ = new QGroupBox(Point_set_base);
        table_frame_->setObjectName(QString::fromUtf8("table_frame_"));

        vboxLayout->addWidget(table_frame_);

        Spacer10_2 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout->addItem(Spacer10_2);

        standardize_options_box_ = new QGroupBox(Point_set_base);
        standardize_options_box_->setObjectName(QString::fromUtf8("standardize_options_box_"));
        vboxLayout1 = new QVBoxLayout(standardize_options_box_);
        vboxLayout1->setSpacing(0);
        vboxLayout1->setContentsMargins(0, 0, 0, 0);
        vboxLayout1->setObjectName(QString::fromUtf8("vboxLayout1"));
        hboxLayout5 = new QHBoxLayout();
        hboxLayout5->setSpacing(6);
        hboxLayout5->setContentsMargins(0, 0, 0, 0);
        hboxLayout5->setObjectName(QString::fromUtf8("hboxLayout5"));
        textLabel1_3 = new QLabel(standardize_options_box_);
        textLabel1_3->setObjectName(QString::fromUtf8("textLabel1_3"));
        sizePolicy.setHeightForWidth(textLabel1_3->sizePolicy().hasHeightForWidth());
        textLabel1_3->setSizePolicy(sizePolicy);
        textLabel1_3->setFont(font);
        textLabel1_3->setWordWrap(false);

        hboxLayout5->addWidget(textLabel1_3);

        line1_3 = new QFrame(standardize_options_box_);
        line1_3->setObjectName(QString::fromUtf8("line1_3"));
        sizePolicy1.setHeightForWidth(line1_3->sizePolicy().hasHeightForWidth());
        line1_3->setSizePolicy(sizePolicy1);
        line1_3->setFrameShape(QFrame::HLine);
        line1_3->setFrameShadow(QFrame::Raised);

        hboxLayout5->addWidget(line1_3);


        vboxLayout1->addLayout(hboxLayout5);

        spacer10 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        vboxLayout1->addItem(spacer10);

        hboxLayout6 = new QHBoxLayout();
        hboxLayout6->setSpacing(6);
        hboxLayout6->setContentsMargins(0, 0, 0, 0);
        hboxLayout6->setObjectName(QString::fromUtf8("hboxLayout6"));
        spacer8_2 = new QSpacerItem(51, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout6->addItem(spacer8_2);

        standardize_checkbox_ = new QCheckBox(standardize_options_box_);
        standardize_checkbox_->setObjectName(QString::fromUtf8("standardize_checkbox_"));

        hboxLayout6->addWidget(standardize_checkbox_);

        spacer8 = new QSpacerItem(51, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout6->addItem(spacer8);


        vboxLayout1->addLayout(hboxLayout6);


        vboxLayout->addWidget(standardize_options_box_);

        spacer5 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);

        vboxLayout->addItem(spacer5);

        QWidget::setTabOrder(num_lags_, lag_sep_);
        QWidget::setTabOrder(lag_sep_, lag_tol_);
        QWidget::setTabOrder(lag_tol_, directions_count_);
        QWidget::setTabOrder(directions_count_, standardize_checkbox_);
        QWidget::setTabOrder(standardize_checkbox_, load_button_);
        QWidget::setTabOrder(load_button_, save_button_);

        retranslateUi(Point_set_base);

        QMetaObject::connectSlotsByName(Point_set_base);
    } // setupUi

    void retranslateUi(QWidget *Point_set_base)
    {
        Point_set_base->setWindowTitle(QApplication::translate("Point_set_base", "Point_set_base", 0, QApplication::UnicodeUTF8));
        load_button_->setText(QApplication::translate("Point_set_base", "Load Parameters...", 0, QApplication::UnicodeUTF8));
        save_button_->setText(QApplication::translate("Point_set_base", "Save", 0, QApplication::UnicodeUTF8));
        textLabel1->setText(QApplication::translate("Point_set_base", "Lags", 0, QApplication::UnicodeUTF8));
        TextLabel4->setText(QApplication::translate("Point_set_base", "Lag tolerance", 0, QApplication::UnicodeUTF8));
        TextLabel1_2->setText(QApplication::translate("Point_set_base", "Lag separation", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("Point_set_base", "Number of lags", 0, QApplication::UnicodeUTF8));
        textLabel1_2->setText(QApplication::translate("Point_set_base", "Directions", 0, QApplication::UnicodeUTF8));
        textLabel1_4->setText(QApplication::translate("Point_set_base", "Number of directions", 0, QApplication::UnicodeUTF8));
        textLabel1_5->setText(QApplication::translate("Point_set_base", "<i>Angles are in degrees.<br>\n"
"Use a tolerance tol > 90 to indicate an omni-directional variogram</i>", 0, QApplication::UnicodeUTF8));
        table_frame_->setTitle(QString());
        standardize_options_box_->setTitle(QString());
        textLabel1_3->setText(QApplication::translate("Point_set_base", "Misc.", 0, QApplication::UnicodeUTF8));
        standardize_checkbox_->setText(QApplication::translate("Point_set_base", "Standardize by Cov(head, tail)", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Point_set_base: public Ui_Point_set_base {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINT_SET_BASE_H
