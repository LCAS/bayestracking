^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package bayes_tracking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.0.2 (2014-09-03)
------------------
* Adding boost as a build an run dependency.
* Changing Licence to GPL
* Contributors: Christian Dondrup

1.3.0 (2018-09-30)
------------------
* fixed tag check for LABELED (`#24 <https://github.com/LCAS/bayestracking/issues/24>`_)
* Contributors: Marc Hanheide

1.2.0 (2018-09-28)
------------------
* Flexible id associations (`#23 <https://github.com/LCAS/bayestracking/issues/23>`_)
  * fix to use local headers first
  * proper support for tagged tracking
  * added tags tracking example
  * fix bug to only assign if not named yet
* Contributors: Marc Hanheide

1.1.0 (2018-09-04)
------------------
* support for Polar Models (`#22 <https://github.com/LCAS/bayestracking/issues/22>`_)
  * Added 2D Polar Observation Model
  * no default parameter values in redeclared functions
  * proper building
  * Changes for backward-compatibility
  * Backward-compatibility of example code
* Fixed condition logic that caused label checks to always be used (`#21 <https://github.com/LCAS/bayestracking/issues/21>`_)
* Update snapcraft.yaml
* Update snapcraft.yaml
* Update snapcraft.yaml
* Contributors: Marc Hanheide, SeanTasker, Strands JENKINS Daemon

1.0.8 (2017-06-13)
------------------
* added cv_bridge to catkin build (`#18 <https://github.com/LCAS/bayestracking/issues/18>`_)
* Tracking NN including labels (`#13 <https://github.com/LCAS/bayestracking/issues/13>`_)
  * Tracking NN including labels
  without auto formatting
  * Implimented NN_LABELED
  Restored default behaviour implimented another association algorithm
  * Included NN_LABEL in Association Matrix
  * correct mistyped word
  * NNJPDA_LABEL implimented
  * Added seq_size and seq_time
* Contributors: Marc Hanheide, Peter Lightbody

1.0.7 (2016-07-02)
------------------
* Revert "Update on bayestracking" (`#17 <https://github.com/LCAS/bayestracking/issues/17>`_)
* Contributors: Marc Hanheide

1.0.6 (2016-07-02)
------------------
* Update on bayestracking (`#16 <https://github.com/LCAS/bayestracking/issues/16>`_)
  * Added option to specify minimum number of observations, minimum interval between observations for new track creation and a variance limit for track deletion
  * Update on mullitracker.h
  Change the naming "varLimit" -> "stdLimit"
  * Update simple_2d_tracking.cpp
  Change naming "varLimit" -> "stdLimit"
* Merge pull request `#11 <https://github.com/LCAS/bayestracking/issues/11>`_ from pet1330/master
  Adapted to allow 3D tracking along side 2D
* removed 3D from simple example
* removed X_SIZE reference
* re-implemented 2D functionality to allow merge
* 3D version
* Contributors: Christian Dondrup, Peter Lightbody, scosar

1.0.5 (2014-10-29)
------------------
* Merge pull request `#9 <https://github.com/LCAS/bayestracking/issues/9>`_ from LCAS/indigo-fix
  Trusty/Indigo fix
* Trusty/Indigo fix
* Contributors: Christian Dondrup

1.0.4 (2014-10-29)
------------------
* Merge pull request `#8 <https://github.com/LCAS/bayestracking/issues/8>`_ from LCAS/opencv
  Changing dependency from opencv2 to cv_bridge for indigo release.
* Changing dependency from opencv2 to cv_bridge for indigo release.
* Contributors: Christian Dondrup

1.0.3 (2014-09-09)
------------------
* "Fixing" the missing JPDA
  The option of using only JPDA as an association algorithm was a relic from past development and is not supported anymore.
  This change just comments the JPDA in the enum, preventing people from trying to use it.
* Using the correct licence now.
* 1.0.2
* Updated changelog for 1.0.2
* Adding boost as a build an run dependency.
* 1.0.1
* Reverting the version number increase.
* Merge pull request `#5 <https://github.com/LCAS/bayestracking/issues/5>`_ from cdondrup/master
  Release preparations
* Contributors: Christian Dondrup

1.0.1 (2014-09-02)
------------------
* Merge pull request `#4 <https://github.com/cdondrup/bayestracking/issues/4>`_ from cdondrup/master
  Renamed `bayestracking` to `bayes_tracking` to comply with ros naming conventions.
* Renamed `bayestracking` to `bayes_tracking` to comply with ros naming conventions.
* Merge pull request `#2 <https://github.com/cdondrup/bayestracking/issues/2>`_ from cdondrup/master
  Adding optional catkin build
* * Enabling normal cmake build if catkin is not available.
  * Updated package.xml and create install targets for ROS release
  * Refactored README to markdown
* Refactored bayestracking library to be a catkin package
  Created a cleaner file structure by seperating the src from the include files.
* Update README
* Update README
* Now generates cmake config files to enable usage of find_package
* Added the creation of a pkg-config file to make the library easier to use.
* PFilter not working right now in simple_tracking.cpp, but all other tested filters (UKF and EKF) worked very well
* fixed size_t uint issue
* compiles now (but simple_tracking example still gives errors)
* initial import form package
* Contributors: Christian Dondrup, Marc Hanheide
