import logging
import os
import json
import time
from rdkit import Chem

REPLACE_LETTER = {"(": "_", ")": "_", "'": "_"}