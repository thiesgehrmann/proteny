
def makeConf(chromosomes, karyotypeFiles, genesFile, linksFile):

  return TXT_circos_conf % { "chrUnits" : "1",
                             "chromosomes": chromosomes,
                             "karyotype" : ','.join(karyotypeFiles),
                             "gene": genesFile,
                             "link": linksFile }

#edef

############################################################################### 

TXT_circos_conf = """
# circos conf
# Automatically generated by proteny tool
# https://github.com/thiesgehrmann/proteny
karyotype = %(karyotype)s
chromosomes_units           = %(chrUnits)s
chromosomes_display_default = no
chromosomes                 = %(chromosomes)s

# Define the ideogram
<ideogram>
  <spacing>
    default = 0.01r
  </spacing>
  radius    = 0.90r
  thickness = 20p
  fill      = yes
  show_label       = yes
  label_font       = default
  label_radius     = dims(image,radius) - 90p
  label_size       = 30
  label_rotate     = no
  label_parallel   = yes
</ideogram>


# Add the gene regions
<plots>
  <plot>
    file = %(gene)s
    type        = tile
    r1          = 1r
    r0          = 0.92r
    orientation = out

    layers      = 1
    margin      = 0.02u
    thickness   = 20
    padding     = 0

    layers_overflow=collapse
    stroke_thickness = 1
    <rules>
      <rule>
        condition    = var(strand) eq "n"
        stroke_color = red
        color        = red
      </rule>
      <rule>
        condition    = var(strand) eq "p"
        stroke_color = blue
        color        = blue
      </rule>
    </rules>
  </plot>
</plots>


# Add the links
<links>
  <link>
    file = %(link)s
    radius        = 0.92r
    bezier_radius = 0r
    thickness     = 2
    ribbon        = yes
    stroke_color     = vdgrey_a4
    stroke_thickness = 1
    
    <rules>
      flow = continue
      <rule>
        condition = 1
        color     = eval(sprintf("%%d,%%d,%%d,0.5", var(r_int), var(g_int), var(b_int)))
        z         = eval( 99999999999999999 - min(var(SIZE1),var(SIZE2)) )
        #flat     = yes
      </rule>
    </rules>
  </link>
</links>

################################################################
# The remaining content is standard and required. It is imported 
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.
<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>
# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>
<<include etc/colors.conf>>
# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
"""

###############################################################################