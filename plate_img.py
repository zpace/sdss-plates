import numpy as np
import matplotlib.pyplot as plt
import montage_wrapper as m
import aplpy
from gz2tools import download_sloan_im
import astropy.table as table
import mpld3
from mpld3 import plugins, utils
import pandas as pd
import matplotlib as mpl


def retrieve_ugr(ra, dec, plate_num, width=1., scale=1.):
    '''
    use montage to get a fits header, and then 3-band fits images
    '''

    ra, dec = str(ra), str(dec)
    plate_num = str(plate_num)
    loc_sky = '{0} {1}'.format(ra, dec)

    m.mHdr(object_or_location=loc_sky, system='equatorial',
           width=width, pix_size=scale, out_file='temp.hdr')
    m.mExec(survey='SDSS', band='u', keep=False, corners=True,
            output_image=plate_num + '_u.fits', region_header='temp.hdr')
    m.mExec(survey='SDSS', band='g', keep=False, corners=True,
            output_image=plate_num + '_g.fits', region_header='temp.hdr')
    m.mExec(survey='SDSS', band='r', keep=False, corners=True,
            output_image=plate_num + '_r.fits', region_header='temp.hdr')


def make_SDSS_FITSFigure(ra, dec, plate_num, width):
    '''
    aspect ratio is weird with this method, so trying something different
    '''

    ra, dec = str(ra), str(dec)
    plate_num = str(plate_num)

    aplpy.make_rgb_cube([plate_num + '_r.fits', plate_num +
                         '_g.fits', plate_num + '_u.fits'],
                        plate_num + '_cube.fits')
    aplpy.make_rgb_image(plate_num + '_cube.fits',
                         plate_num + '_rgb.png')
    plt.close('all')
    f = aplpy.FITSFigure(plate_num + '_cube_2d.fits', figsize=(10, 10))
    f.show_rgb(plate_num + '_rgb.png')
    f.add_grid()
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='x-small')
    f.axis_labels.set_xtext('Right Ascension (J2000)')
    f.axis_labels.set_ytext('Declination (J2000)')

    plt.tight_layout()
    plt.show()


def make_SDSS_field_cutout(ra, dec, plate_num, width=3.1, scale=5.5):
    '''
    use gz2tools's image downloader

    width argument is *in degrees*, and download_sloan_im takes width and height in pixels (with scale in arcsec)

    NOTE: there seems to be an artificial limit of 2400*2400 pix. Current defaults seem to play nice, but be careful!

    NOTE: this will likely break along the RA=0h meridian! May hack a fix, or transition to basemap at a later date!
    '''

    plt.close('all')
    img = download_sloan_im(ra, dec, scale=scale, width=(width*3600.)//scale,
                            height=(width*3600.)//scale)
    print img.shape, (width*3600.)//scale
    ralims = [sum(i) for i in zip([ra, ] * 2, [width/2., -width/2.])]
    declims = [sum(i) for i in zip([dec, ] * 2, [width/2., -width/2.])]

    fig = plt.figure(figsize=(8.5, 8.), dpi=300)
    ax = fig.add_subplot(111)

    ax.imshow(img, origin='lower', aspect='equal',
              extent=ralims + declims)

    # add in plate shape and center post
    plate = plt.Circle((ra, dec), 1.49, color='w', fill=False)
    ax.add_artist(plate)

    # read in locations of objects on plate
    plate_objs = table.Table.read(str(plate_num) + '.csv', format='csv',
                                  comment='\#')
    stars = plate_objs[plate_objs['class'] == 'STAR']
    QSOs = plate_objs[plate_objs['class'] == 'QSO']
    galaxies = plate_objs[plate_objs['class'] == 'GALAXY']

    # add locations of objects on plate
    ax.scatter(stars['ra'], stars['dec'],
               marker='D', facecolor='none', edgecolor='b',
               label='Stars')
    ax.scatter(QSOs['ra'], QSOs['dec'],
               marker='o', facecolor='none', edgecolor='g',
               label='QSOs')
    ax.scatter(galaxies['ra'], galaxies['dec'],
               marker='+', facecolor='none', edgecolor='r',
               label='Galaxies')
    plt.legend(loc='best')

    #ax.legend(loc='upper right')

    # make everything look nice
    ax.set_xlabel('RA [deg]', size=20)
    ax.set_ylabel('Dec [deg]', size=20)
    ax.set_title('Plate {}'.format(plate_num), size=20)
    ax.set_xlim(ralims)
    ax.set_ylim(declims)

    ax.tick_params(axis='both', colors='white')
    ax.xaxis.set_tick_params(width=4, length=12)
    ax.yaxis.set_tick_params(width=4, length=12)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_color('k')
        tick.label.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_color('k')
        tick.label.set_fontsize(20)

    plt.tight_layout()
    plt.savefig(str(plate_num) + '.png')


class ClickInfo(mpld3.plugins.PluginBase):

    """mpld3 Plugin for getting info on click        """

    JAVASCRIPT = """
    mpld3.register_plugin("clickinfo", ClickInfo);
    ClickInfo.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo.prototype.constructor = ClickInfo;
    ClickInfo.prototype.requiredProps = ["id", "urls"];
    function ClickInfo(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    ClickInfo.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        urls = this.props.urls;
        obj.elements().on("mousedown",
                          function(d, i){
                            window.open(urls[i], '_blank')});
    }
    """

    def __init__(self, points, urls):
        self.points = points
        self.urls = urls
        if isinstance(points, mpl.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None
        self.dict_ = {"type": "clickinfo",
                      "id": utils.get_id(points, suffix),
                      "urls": urls}


def make_SDSS_field_mpld3(ra, dec, plate_num, width=3.1, scale=5.5):
    '''
    use gz2tools's image downloader to show SDSS field with fibers overlaid

    on click, opens the corresponding SDSS Explore page

    width argument is *in degrees*, and download_sloan_im takes width and height in pixels (with scale in arcsec)

    NOTE: there seems to be an artificial limit of 2400*2400 pix. Current defaults seem to play nice, but be careful!

    NOTE: this will likely break along the RA=0h meridian! May hack a fix, or transition to basemap at a later date!
    '''

    # Define some CSS to control our custom labels
    css = """
    table
    {
      border-collapse: collapse;
    }
    th
    {
      color: #ffffff;
      background-color: #000000;
    }
    td
    {
      background-color: #cccccc;
    }
    table, th, td
    {
      font-family:Arial, Helvetica, sans-serif;
      border: 1px solid black;
      text-align: right;
    }
    """

    plt.close('all')
    img = download_sloan_im(ra, dec, scale=scale, width=(width*3600.)//scale,
                            height=(width*3600.)//scale)
    print img.shape, (width*3600.)//scale
    ralims = [sum(i) for i in zip([ra, ] * 2, [width/2., -width/2.])]
    declims = [sum(i) for i in zip([dec, ] * 2, [width/2., -width/2.])]

    fig = plt.figure(figsize=(8.5, 8.), dpi=100)
    ax = fig.add_subplot(111)

    ax.imshow(img, origin='upper', aspect='equal',
              extent=ralims + declims)

    # add in plate shape and center post
    plate = plt.Circle((ra, dec), 1.49, color='w', fill=False)
    ax.add_artist(plate)

    # read in locations of objects on plate
    plate_objs = table.Table.read(str(plate_num) + '.csv', format='csv',
                                  comment='\#')
    plate_objs = plate_objs[plate_objs['targetObjID'] != 0]
    stars = plate_objs[plate_objs['class'] == 'STAR']
    QSOs = plate_objs[plate_objs['class'] == 'QSO']
    galaxies = plate_objs[plate_objs['class'] == 'GALAXY']

    # add locations of objects on plate

    stars_points = ax.scatter(stars['ra'], stars['dec'],
                              marker='D', edgecolors='b', facecolors='none',
                              label='Stars')
    QSOs_points = ax.scatter(QSOs['ra'], QSOs['dec'],
                             marker='o', edgecolors='g', facecolors='none',
                             label='QSOs')
    galaxies_points = ax.scatter(galaxies['ra'], galaxies['dec'],
                                 marker='v', edgecolors='r', facecolors='none',
                                 label='Galaxies')

    # ax.legend(loc='best') #not implemented in mpld3 yet
    # replaced with figtext
    xtext =
    ytext =
    ax.figtext()

    # make everything look nice
    ax.set_xlabel('RA [deg]', size=20)
    ax.set_ylabel('Dec [deg]', size=20)
    ax.set_title('Plate {}'.format(plate_num), size=20)
    ax.set_xlim(ralims)
    ax.set_ylim(declims)

    ax.tick_params(axis='both', colors='white')
    ax.xaxis.set_tick_params(width=4, length=12)
    ax.yaxis.set_tick_params(width=4, length=12)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_color('k')
        tick.label.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_color('k')
        tick.label.set_fontsize(20)

    plt.tight_layout()

    # read in everything as dataframes, since that makes labeling easier
    plate_objs_df = pd.read_csv(str(plate_num) + '.csv', comment='#')
    plate_objs_df = plate_objs_df[plate_objs_df['targetObjID'] != 0]

    stars_df = plate_objs_df[plate_objs_df['class'] == 'STAR']
    QSOs_df = plate_objs_df[plate_objs_df['class'] == 'QSO']
    galaxies_df = plate_objs_df[plate_objs_df['class'] == 'GALAXY']

    labels_stars = []
    for i in range(len(stars)):
        label = stars_df.iloc[i].T
        label = pd.DataFrame({'Row {0}'.format(i): label})
        # .to_html() is unicode; so make leading 'u' go away with str()
        labels_stars.append(str(label.to_html()))

    labels_QSOs = []
    for i in range(len(QSOs)):
        label = QSOs_df.iloc[i].T
        label = pd.DataFrame({'Row {0}'.format(i): label})
        # .to_html() is unicode; so make leading 'u' go away with str()
        labels_QSOs.append(str(label.to_html()))

    labels_galaxies = []
    for i in range(len(galaxies)):
        label = galaxies_df.iloc[i].T
        label = pd.DataFrame({'Row {0}'.format(i): label})
        # .to_html() is unicode; so make leading 'u' go away with str()
        labels_galaxies.append(str(label.to_html()))

    url_base = \
        'http://skyserver.sdss.org/dr12/en/tools/explore/summary.aspx?id='
    urls_stars = [url_base + str(i) for i in stars['targetObjID']]
    urls_QSOs = [url_base + str(i) for i in QSOs['targetObjID']]
    urls_galaxies = [url_base + str(i) for i in galaxies['targetObjID']]

    # instantiate the plugins, one for each dataset
    tooltip_stars = plugins.PointHTMLTooltip(stars_points, labels_stars,
                                             voffset=10, hoffset=10, css=css)
    tooltip_QSOs = plugins.PointHTMLTooltip(QSOs_points, labels_QSOs,
                                            voffset=10, hoffset=10, css=css)
    tooltip_galaxies = plugins.PointHTMLTooltip(
        galaxies_points, labels_galaxies, voffset=10, hoffset=10, css=css)
    # link everything
    plugins.connect(fig, tooltip_stars)
    plugins.connect(fig, tooltip_QSOs)
    plugins.connect(fig, tooltip_galaxies)

    plugins.connect(fig, ClickInfo(stars_points, urls_stars))
    plugins.connect(fig, ClickInfo(QSOs_points, urls_QSOs))
    plugins.connect(fig, ClickInfo(galaxies_points, urls_galaxies))

    mpld3.save_html(fig, str(plate_num) + '.html')
