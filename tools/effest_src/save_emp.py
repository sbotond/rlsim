#
# Copyright (C) 2013 EMBL - European Bioinformatics Institute
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# Neither the institution name nor the name rlsim
# can be used to endorse or promote products derived from
# this software without prior written permission. For
# written permission, please contact <sbotond@ebi.ac.uk>.

# Products derived from this software may not be called
# rlsim nor may rlsim appear in their
# names without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

class SaveEmp:
    def __init__(self, fname):
        self.fname  = fname
        self.fh     = open(fname, "w")
        self.pool   = OrderedDict() 

    def add_obj(self, name, obj):
        self.pool[name] = obj

    def add_kv(self, name, k, v):
        if len(k) != len(v):
            L.Fatal("Number of keys is not equalt to the number of values!")
        tmp = OrderedDict() 
        for i in xrange(len(k)):
           tmp[str(k[i])]   = v[i] 
        self.add_obj(name, tmp)

    def save_objs(self):
        data    = json.dumps(self.pool, indent=True)
        self.fh.write(data)
        self.fh.flush()
        self.fh.close()
