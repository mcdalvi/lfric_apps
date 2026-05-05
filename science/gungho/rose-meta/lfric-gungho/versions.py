import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""

class vn3.1_t474(MacroUpgrade):
    # Upgrade macro for #474 by Mohit Dalvi
    # Add namelist and related items for Nudging functionality

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t474"

    def upgrade(self, config, meta_config=None):
        # Add settings - new items in files and multires_coupling
        self.add_setting( config,
            ["namelist:files", "nudging_directory"],"''")
        self.add_setting( config,
            ["namelist:files", "nudging_filename"],"''")

        self.add_setting( config,
            ["namelist:multires_coupling", "coarse_nudging"], ".false," )
        self.add_setting( config,
            ["namelist:multires_coupling", "nudging_mesh_name"], "''")

        # Add new nudging namelist 
        # Append after 'multires_coupling' in configuration.nml
        source = self.get_setting_value( config,
                ["file:configuration.nml","source"] )
        source = re.sub(r'namelist:multires_coupling',
                 r'namelist:multires_coupling' + '\n' + ' (namelist:nudging)',
                 source)
        self.change_setting_value( config,
                ["file:configuration.nml","source"], source )

        self.add_setting(config, ["namelist:nudging"])
        self.add_setting(config, ["namelist:nudging","nudge_data_levels"], 0 )
        self.add_setting(config,
                         ["namelist:nudging","nudging_level_bottom"], 0 )
        self.add_setting(config,
                         ["namelist:nudging","nudging_level_top"], 0 )
        self.add_setting(config, ["namelist:nudging","nudge_source"], "''" )
        self.add_setting(config,
                         ["namelist:nudging","nudging_width_bottom"], 0 )
        self.add_setting(config,
                         ["namelist:nudging","nudging_width_top"], 0 )

        return config, self.reports
