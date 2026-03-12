#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Interactive RSS image and spectrum explorer.

Displays a spatial flux image from an RSS FITS file (lvmSFrame,
lvmCFrame, or rss_combine output) with two interactive sliders
and a click-to-spectrum panel.  Can be used from a Jupyter notebook
or run directly from the command line.

Command line usage:

    rss_explore.py filename [--center WAVELENGTH] [--width WIDTH] [--ext EXT]

    filename
        RSS FITS file (lvmSFrame, lvmCFrame, or rss_combine output).

    --center WAVELENGTH
        Starting center wavelength in Angstroms (default: 6563, Ha).

    --width WIDTH
        Starting wavelength window width in Angstroms (default: 20).

    --ext EXT
        FITS extension to display (default: FLUX).

    Examples:

        rss_explore.py lvmSFrame-00007373.fits
        rss_explore.py lvmSFrame-00007373.fits --center 5007 --width 15
        rss_explore.py lvmSFrame-00007373.fits --ext IVAR

Usage (in a Jupyter notebook):

    %matplotlib widget          # must be the first line, before imports

    from rss_explore import RSSExplorer

    exp = RSSExplorer('myfile.fits')
    exp.show()                       # defaults: center=6563 A, width=20 A
    exp.show(center=5007, width=15)  # start on [O III]

    Multiple independent calls to show() can be active simultaneously;
    each notebook cell has its own sliders, image, and spectrum panel.

Interactions (both notebook and command-line):

    Sliders
        Center (A)  - center wavelength of the displayed passband.
        Width  (A)  - width of the passband in Angstroms.
        Moving either slider recomputes and redraws the spatial image.

    Click on the image
        Selects the nearest fibre.  Its full spectrum is displayed in
        the lower panel with a red cross marking the selected position.
        The panel title shows the fibre number, RA, and Dec in degrees.

    Drag a box on the spectrum
        Zooms the spectrum panel to the selected wavelength and flux
        region.  When zoomed, switching to a different fibre keeps the
        same wavelength and flux range so fibres can be compared
        directly on the same scale.

    Z button (command line) / Reset zoom button (notebook)
        Restores the spectrum panel to the full wavelength range.

Description:

The RSSExplorer class loads the FITS file and precomputes the WCS,
per-fibre pixel positions, and per-fibre pixel masks once on
initialisation.  This makes the one-time startup take a few seconds
but keeps all subsequent slider and click interactions fast (< 1 s).

Each call to show() creates its own figure and event handlers as
local closures so that multiple show() calls in a notebook are fully
independent.

The command-line version (show_cli) uses matplotlib.widgets.Slider
and matplotlib.widgets.Button in a standalone figure window.  The
notebook version (show) uses ipywidgets sliders embedded alongside
a %matplotlib widget canvas.

Requires:
    ipywidgets  - for notebook sliders
    ipympl      - for interactive notebook canvas (%matplotlib widget)

    pip install ipywidgets ipympl

The routine reuses helper functions from rss2image.py
(get_flux, clean_slitmap, create_wcs_from_points,
create_wcs_manually, get_redshift_correction).

Notes:

The RA/Dec tick labels on the image are approximate: RA values are
computed at the middle row of the image and Dec values at the middle
column.  The error is negligible over a typical LVM field (~0.5 deg).

History:

250220 ksl Coding begun, based on rss2image.py.

'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector, Button as MplButton
import astropy.io.fits as fits
from astropy.table import Table, join
from ipywidgets import FloatSlider, VBox
from IPython.display import display

from rss2image import (
    get_flux,
    clean_slitmap,
    create_wcs_from_points,
    create_wcs_manually,
    get_redshift_correction,
)


class RSSExplorer:
    '''
    Interactive RSS image explorer.

    Load a FITS RSS file once, then explore the spatial flux
    distribution interactively across wavelength via sliders.
    Click on any fibre in the image to display its spectrum below.

    Parameters
    ----------
    filename : str
        Path to the RSS FITS file.
    xscale : float, optional
        Pixel scale in arcsec/pixel used when building a manual WCS.
        Default 1.
    xsize : float, optional
        Approximate field size in degrees used for the manual WCS.
        Default 0.5.
    ext : str, optional
        FITS extension to read flux from.  Default 'FLUX'.
    '''

    def __init__(self, filename, xscale=1, xsize=0.5, ext='FLUX'):
        self.filename = filename
        self.xscale = xscale
        self.xsize = xsize
        self.ext = ext

        print(f'Loading {filename} ...')
        self.rss = fits.open(filename)

        phead = self.rss['PRIMARY'].header
        xtab = Table(self.rss['SLITMAP'].data)
        xtab['Row'] = np.arange(len(xtab))
        self.qtab = clean_slitmap(xtab)

        # Build WCS - try X/Y columns first, fall back to header RA/Dec
        have_wcs = False
        try:
            xpos = self.qtab['X']
            ypos = self.qtab['Y']
            wcs = create_wcs_from_points(xpos, ypos, self.qtab['ra'], self.qtab['dec'])
            have_wcs = True
        except Exception:
            pass

        if not have_wcs:
            try:
                xra = phead['POSCIRA']
                xdec = phead['POSCIDE']
                wcs = create_wcs_manually(xra, xdec, xscale, xsize)
                have_wcs = True
                print('Generated WCS from primary header')
            except Exception:
                pass

        if not have_wcs:
            raise RuntimeError(f'Could not generate WCS from {filename}')

        self.wcs = wcs

        xra = wcs.wcs.crval[0]
        xdec = wcs.wcs.crval[1]
        self.redshift_correction = get_redshift_correction(xra, xdec)

        # Precompute pixel positions for every fibre
        x_pix, y_pix = wcs.world_to_pixel_values(self.qtab['ra'], self.qtab['dec'])
        self.qtab['X'] = x_pix
        self.qtab['Y'] = y_pix

        npix = int(np.nanmax([np.nanmax(x_pix), np.nanmax(y_pix)])) + 1
        self.npix = npix
        rspaxel = 35.0 / 2.0 / xscale   # fibre radius in pixels

        # Precompute pixel masks for every fibre (the slow loop, done once).
        # _fiber_masks maps Row -> (row_indices, col_indices) into the image.
        print('Precomputing fibre pixel masks ...')
        xima, yima = np.meshgrid(np.arange(npix), np.arange(npix))
        self._fiber_masks = {}
        for row in self.qtab:
            sel = (xima - row['X'])**2 + (yima - row['Y'])**2 <= rspaxel**2
            self._fiber_masks[int(row['Row'])] = np.where(sel)

        # Wavelength array and limits
        self._wave = self.rss['WAVE'].data
        self.wave_min = float(np.nanmin(self._wave))
        self.wave_max = float(np.nanmax(self._wave))

        print(f'Ready.  Wavelength range: {self.wave_min:.1f} - {self.wave_max:.1f} A')

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_image(self, center, width):
        '''
        Extract flux in [center-width/2, center+width/2] and paint
        the spatial image using precomputed fibre masks.

        Returns
        -------
        ima : 2-D ndarray
        wmin, wmax : float  (redshift-corrected wavelength limits used)
        '''
        rc = self.redshift_correction
        wmin = (center - width / 2.0) * rc
        wmax = (center + width / 2.0) * rc

        if self.ext == 'FLUX':
            _, xflux = get_flux(self.rss, wmin, wmax, xtype='sum', ext=self.ext)
        else:
            _, xflux = get_flux(self.rss, wmin, wmax, xtype='med', ext=self.ext)

        merged = join(self.qtab, xflux, keys='Row', join_type='left')

        ima = np.full((self.npix, self.npix), np.nan)
        for row in merged:
            mask = self._fiber_masks.get(int(row['Row']))
            if mask is not None:
                ima[mask] = row['FLUX']

        return ima, wmin, wmax

    def _set_wcs_ticks(self, ax, n_ticks=5):
        '''
        Label the image axes with approximate RA and Dec values.

        RA ticks are computed along the middle row of the image;
        Dec ticks along the middle column.  Both are approximate
        because the projection is not strictly separable, but the
        error is negligible over a small LVM field.
        '''
        npix = self.npix
        mid  = npix / 2.0
        tick_pix = np.linspace(0, npix - 1, n_ticks)

        ra_vals, _  = self.wcs.pixel_to_world_values(tick_pix, np.full(n_ticks, mid))
        _, dec_vals = self.wcs.pixel_to_world_values(np.full(n_ticks, mid), tick_pix)

        ax.set_xticks(tick_pix)
        ax.set_xticklabels([f'{r:.2f}' for r in ra_vals])
        ax.set_yticks(tick_pix)
        ax.set_yticklabels([f'{d:.2f}' for d in dec_vals])

    @staticmethod
    def _get_clim(ima):
        finite = ima[np.isfinite(ima)]
        if len(finite) > 0:
            return np.nanpercentile(finite, [5, 95])
        return 0.0, 1.0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def show(self, center=6563.0, width=20.0):
        '''
        Launch the interactive explorer in a Jupyter notebook cell.

        Each call to show() is fully independent: its own figure,
        sliders, and event handlers do not interact with other calls.

        Requires %matplotlib widget (ipympl) to be active.

        Parameters
        ----------
        center : float, optional
            Starting center wavelength in Angstroms.  Default 6563 (Ha).
        width : float, optional
            Starting wavelength window width in Angstroms.  Default 20.
        '''
        wlo = self.wave_min
        whi = self.wave_max

        center = float(np.clip(center, wlo, whi))
        width  = float(np.clip(width,  1.0, min(200.0, whi - wlo)))

        # Create the figure with ioff() so it does not auto-display;
        # we embed it explicitly in the VBox below.
        # Extra bottom margin leaves room for the Reset zoom button.
        with plt.ioff():
            fig, (ax_img, ax_spec) = plt.subplots(
                2, 1, figsize=(8, 10),
                gridspec_kw={'height_ratios': [1.4, 1]},
            )
        fig.subplots_adjust(bottom=0.15)
        fig.canvas.toolbar_visible = False

        # Initial image
        ima, wmin, wmax = self._make_image(center, width)
        vmin, vmax = self._get_clim(ima)
        imshow = ax_img.imshow(ima, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        plt.colorbar(imshow, ax=ax_img, label=self.ext)
        ax_img.set_title(
            f'{self.ext}   {wmin:.1f} - {wmax:.1f} A'
            f'  (center {center:.1f}, width {width:.1f} A)'
        )
        ax_img.set_xlabel('RA')
        ax_img.set_ylabel('Dec')
        self._set_wcs_ticks(ax_img)

        # Marker for selected fibre (hidden until first click)
        marker, = ax_img.plot([], [], 'r+', markersize=18, markeredgewidth=2, zorder=5)

        # Spectrum panel (empty until first click)
        spec_line, = ax_spec.plot([], [], 'b-', linewidth=0.8)
        ax_spec.set_xlabel('Wavelength (A)')
        ax_spec.set_ylabel(self.ext)
        ax_spec.set_title('Click on the image to show a spectrum')

        fig.tight_layout()

        # ----------------------------------------------------------
        # Event handlers: local closures over this call's local vars.
        # This ensures multiple show() calls are fully independent.
        # ----------------------------------------------------------

        # State shared across the closures below.
        # Use single-element lists so closures can rebind the value.
        current_spec = [None]   # most recently displayed spectrum array
        is_zoomed    = [False]  # True once user has drawn a zoom box

        def _set_full_limits(spec):
            '''Show the full wavelength range, y scaled to the spectrum.'''
            valid = np.isfinite(spec)
            if not np.any(valid):
                return
            ymin = float(np.nanmin(spec[valid]))
            ymax = float(np.nanmax(spec[valid]))
            margin = 0.05 * (ymax - ymin) if ymax != ymin else max(abs(ymax) * 0.05, 1e-30)
            ax_spec.set_xlim(self._wave[0], self._wave[-1])
            ax_spec.set_ylim(ymin - margin, ymax + margin)

        def on_click(event):
            # Image axes: single click selects a fibre
            if event.inaxes is ax_img and not event.dblclick:
                if event.xdata is None or event.ydata is None:
                    return

                dist2 = (self.qtab['X'] - event.xdata)**2 + (self.qtab['Y'] - event.ydata)**2
                idx       = int(np.argmin(dist2))
                row_num   = int(self.qtab['Row'][idx])
                fiber_x   = float(self.qtab['X'][idx])
                fiber_y   = float(self.qtab['Y'][idx])
                fiber_ra  = float(self.qtab['ra'][idx])
                fiber_dec = float(self.qtab['dec'][idx])

                marker.set_data([fiber_x], [fiber_y])

                spec = self.rss[self.ext].data[row_num, :]
                current_spec[0] = spec
                spec_line.set_data(self._wave, spec)

                if not is_zoomed[0]:
                    _set_full_limits(spec)
                # When zoomed, leave both x and y limits untouched so
                # successive fibres are directly comparable on the same scale.

                ax_spec.set_title(
                    f'Fiber {row_num}    RA {fiber_ra:.4f}    Dec {fiber_dec:.4f}'
                )
                fig.canvas.draw_idle()

        def on_box_select(eclick, erelease):
            '''Zoom the spectrum axes to the drawn box.'''
            x1, x2 = sorted([eclick.xdata, erelease.xdata])
            y1, y2 = sorted([eclick.ydata, erelease.ydata])
            if x1 is None or x2 is None or x1 == x2:
                return
            ax_spec.set_xlim(x1, x2)
            ax_spec.set_ylim(y1, y2)
            is_zoomed[0] = True
            fig.canvas.draw_idle()

        def on_slider_change(change):
            c = center_slider.value
            w = width_slider.value
            ima, wmin, wmax = self._make_image(c, w)
            vmin, vmax = self._get_clim(ima)
            imshow.set_data(ima)
            imshow.set_clim(vmin, vmax)
            ax_img.set_title(
                f'{self.ext}   {wmin:.1f} - {wmax:.1f} A'
                f'  (center {c:.1f}, width {w:.1f} A)'
            )
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', on_click)

        # RectangleSelector on the spectrum axes for box zoom.
        # Stored on fig to prevent garbage collection.
        fig._selector = RectangleSelector(
            ax_spec, on_box_select,
            useblit=False,        # True snapshots the empty background on init
            button=[1],           # left mouse button only
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=False,
        )

        # Sliders
        center_slider = FloatSlider(
            value=center, min=wlo, max=whi, step=0.5,
            description='Center (A):',
            continuous_update=False,
            layout={'width': '75%'},
            style={'description_width': '100px'},
        )
        width_slider = FloatSlider(
            value=width, min=1.0, max=200.0, step=1.0,
            description='Width (A):',
            continuous_update=False,
            layout={'width': '75%'},
            style={'description_width': '100px'},
        )

        center_slider.observe(on_slider_change, names='value')
        width_slider.observe(on_slider_change, names='value')

        # Reset zoom button: a matplotlib Button placed inside the figure
        # just below the spectrum panel, so it sits close to the plot.
        # Stored on fig to prevent garbage collection.
        ax_btn = fig.add_axes([0.02, 0.02, 0.18, 0.03])
        fig._reset_btn = MplButton(ax_btn, 'Reset zoom')

        def on_reset(event):
            if current_spec[0] is not None:
                is_zoomed[0] = False
                _set_full_limits(current_spec[0])
                fig.canvas.draw_idle()

        fig._reset_btn.on_clicked(on_reset)

        # Embed sliders and canvas together so this cell is self-contained
        # and independent of other show() calls.
        display(VBox([center_slider, width_slider, fig.canvas]))

    # ------------------------------------------------------------------

    def show_cli(self, center=6563.0, width=20.0):
        '''
        Launch the interactive explorer as a standalone matplotlib window
        for use when running from the command line.

        Parameters
        ----------
        center : float, optional
            Starting center wavelength in Angstroms.  Default 6563 (Ha).
        width : float, optional
            Starting wavelength window width in Angstroms.  Default 20.
        '''
        from matplotlib.widgets import Slider

        wlo = self.wave_min
        whi = self.wave_max

        center = float(np.clip(center, wlo, whi))
        width  = float(np.clip(width,  1.0, min(200.0, whi - wlo)))

        fig = plt.figure(figsize=(10, 10))

        # Axes positions [left, bottom, width, height] in figure fractions.
        #
        #  0.96 ┐
        #       │  ax_img  (image)
        #  0.49 ┘
        #         0.08 gap for axis labels / spectrum title
        #  0.41 ┐
        #       │  ax_spec (spectrum)   ax_btn ("Z") at top-right corner
        #  0.22 ┘
        #  0.15 ┐  ax_sc  (center slider, width 0.70)
        #  0.12 ┘
        #  0.09 ┐  ax_sw  (width slider, width 0.70)
        #  0.06 ┘
        ax_img  = fig.add_axes([0.08, 0.49, 0.87, 0.47])
        ax_spec = fig.add_axes([0.08, 0.22, 0.82, 0.19])
        ax_sc   = fig.add_axes([0.165, 0.12, 0.70, 0.03])    # center slider
        ax_sw   = fig.add_axes([0.165, 0.06, 0.70, 0.03])    # width slider
        ax_btn  = fig.add_axes([0.91, 0.37, 0.04, 0.03])    # "Z" reset button

        # Initial image
        ima, wmin, wmax = self._make_image(center, width)
        vmin, vmax = self._get_clim(ima)
        imshow = ax_img.imshow(ima, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        fig.colorbar(imshow, ax=ax_img, label=self.ext)
        ax_img.set_title(
            f'{self.ext}   {wmin:.1f} - {wmax:.1f} A'
            f'  (center {center:.1f}, width {width:.1f} A)'
        )
        ax_img.set_xlabel('RA')
        ax_img.set_ylabel('Dec')
        self._set_wcs_ticks(ax_img)

        marker, = ax_img.plot([], [], 'r+', markersize=18, markeredgewidth=2, zorder=5)

        spec_line, = ax_spec.plot([], [], 'b-', linewidth=0.8)
        ax_spec.set_xlabel('Wavelength (A)')
        ax_spec.set_ylabel(self.ext)
        ax_spec.set_title('Click on the image to show a spectrum')

        s_center = Slider(ax_sc, 'Center (A)', wlo, whi,    valinit=center, valstep=0.5)
        s_width  = Slider(ax_sw, 'Width (A)',  1.0, 200.0, valinit=width,  valstep=1.0)
        reset_btn = MplButton(ax_btn, 'Z')

        # ----------------------------------------------------------
        # State and event handlers (identical logic to show())
        # ----------------------------------------------------------

        current_spec = [None]
        is_zoomed    = [False]

        def _set_full_limits(spec):
            valid = np.isfinite(spec)
            if not np.any(valid):
                return
            ymin = float(np.nanmin(spec[valid]))
            ymax = float(np.nanmax(spec[valid]))
            margin = 0.05 * (ymax - ymin) if ymax != ymin else max(abs(ymax) * 0.05, 1e-30)
            ax_spec.set_xlim(self._wave[0], self._wave[-1])
            ax_spec.set_ylim(ymin - margin, ymax + margin)

        def on_slider_change(val):
            c = s_center.val
            w = s_width.val
            ima, wmin, wmax = self._make_image(c, w)
            vmin, vmax = self._get_clim(ima)
            imshow.set_data(ima)
            imshow.set_clim(vmin, vmax)
            ax_img.set_title(
                f'{self.ext}   {wmin:.1f} - {wmax:.1f} A'
                f'  (center {c:.1f}, width {w:.1f} A)'
            )
            fig.canvas.draw_idle()

        def on_click(event):
            if event.inaxes is ax_img and not event.dblclick:
                if event.xdata is None or event.ydata is None:
                    return
                dist2 = (self.qtab['X'] - event.xdata)**2 + (self.qtab['Y'] - event.ydata)**2
                idx       = int(np.argmin(dist2))
                row_num   = int(self.qtab['Row'][idx])
                fiber_x   = float(self.qtab['X'][idx])
                fiber_y   = float(self.qtab['Y'][idx])
                fiber_ra  = float(self.qtab['ra'][idx])
                fiber_dec = float(self.qtab['dec'][idx])

                marker.set_data([fiber_x], [fiber_y])

                spec = self.rss[self.ext].data[row_num, :]
                current_spec[0] = spec
                spec_line.set_data(self._wave, spec)

                if not is_zoomed[0]:
                    _set_full_limits(spec)

                ax_spec.set_title(
                    f'Fiber {row_num}    RA {fiber_ra:.4f}    Dec {fiber_dec:.4f}'
                )
                fig.canvas.draw_idle()

        def on_box_select(eclick, erelease):
            x1, x2 = sorted([eclick.xdata, erelease.xdata])
            y1, y2 = sorted([eclick.ydata, erelease.ydata])
            if x1 is None or x2 is None or x1 == x2:
                return
            ax_spec.set_xlim(x1, x2)
            ax_spec.set_ylim(y1, y2)
            is_zoomed[0] = True
            fig.canvas.draw_idle()

        def on_reset(event):
            if current_spec[0] is not None:
                is_zoomed[0] = False
                _set_full_limits(current_spec[0])
                fig.canvas.draw_idle()

        s_center.on_changed(on_slider_change)
        s_width.on_changed(on_slider_change)
        fig.canvas.mpl_connect('button_press_event', on_click)
        reset_btn.on_clicked(on_reset)

        fig._selector = RectangleSelector(
            ax_spec, on_box_select,
            useblit=False,
            button=[1],
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=False,
        )
        fig._reset_btn = reset_btn

        plt.show()


# ======================================================================
# Command-line interface
# ======================================================================

def steer(argv):
    '''Parse arguments and launch the CLI explorer.'''
    import argparse

    parser = argparse.ArgumentParser(
        description='Interactive RSS image and spectrum explorer.'
    )
    parser.add_argument('filename',
                        help='RSS FITS file (lvmSFrame, lvmCFrame, or rss_combine output)')
    parser.add_argument('--center', type=float, default=6563.0,
                        help='Starting center wavelength in Angstroms (default: 6563)')
    parser.add_argument('--width',  type=float, default=20.0,
                        help='Starting wavelength window width in Angstroms (default: 20)')
    parser.add_argument('--ext',    default='FLUX',
                        help='FITS extension to display (default: FLUX)')

    args = parser.parse_args(argv[1:])

    exp = RSSExplorer(args.filename, ext=args.ext)
    exp.show_cli(center=args.center, width=args.width)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
