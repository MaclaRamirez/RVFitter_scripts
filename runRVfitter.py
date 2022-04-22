import sys
import argparse
import matplotlib.pyplot as plt
from RVFitter import RVFitter
import os

plt.switch_backend('QT4Agg')


def clicker(event, mode):
    # global clickValue
    if mode == 'vertical':
        if event.inaxes:
            clickValue['result'] = [event.xdata, event.ydata]
            plt.axvline(x=event.xdata, linestyle='-', color='green')
            plt.draw()
    elif mode == 'horizontal':
        if event.inaxes:
            clickValue['result'] = [event.xdata, event.ydata]
            plt.axhline(y=event.ydata, linestyle='-', color='red')
            plt.draw()
    else:
        print("Unknown position: ", mode)
        print("Choose \"vertical\" or \"horizontal\"")
        raise SystemExit


def on_key(event):
    if event.key == 'enter':
        plt.close()


def close_key(event):
    if event.key == 'y':
        skip['value'] = False
        plt.close()
    elif event.key == 'n':
        skip['value'] = True
        plt.close()
    else:
        print('Unknown key')


def choose_line(rvObject, line):
    print('Do you want to fit this line? (y/n)')
    global skip
    fig, ax = plt.subplots()
    rvObject.plot_line(line=line,
                       ax=ax,
                       title_prefix='Do you want to fit this line? (y/n)\n')
    skip = {}
    skip['value'] = None
    fig.canvas.mpl_connect('key_press_event', lambda event: close_key(event))
    plt.show()
    return skip['value']


def get_clickValue(rvObject,
                   line,
                   position,
                   mode,
                   previous_value=None,
                   fig=None,
                   ax=None):
    if fig == None or ax == None:
        fig, ax = plt.subplots()
    else:
        fig = fig
        ax = ax
    global clickValue
    while True:
        if mode == 'vertical':
            line.plot_normed_spectrum(
                ax=ax,
                title_prefix=
                'Pick {position} boundary wavelength for clipping of\n'.format(
                    position=position))
            line.plot_clips(ax=ax)
            if previous_value is not None:
                ax.axvline(previous_value, color='green')
        elif mode == 'horizontal':
            rvObject.plot_line(
                line=line,
                ax=ax,
                title_prefix=
                'Pick {position} boundary wavelength for normalizing of\n '.
                format(position=position))
            if previous_value is not None:
                ax.axhline(previous_value, color='red')
        else:
            print('mode not known, pick horizontal or vertical')
            raise SystemExit

        clickValue = {}
        clickValue['result'] = [None, None]
        fig.canvas.mpl_connect('button_press_event',
                               lambda event: clicker(event, mode=mode))
        fig.canvas.mpl_connect('key_press_event', on_key)
        plt.show()
        if mode == 'vertical':
            if not clickValue['result'][0] is None:
                foundValue = clickValue['result'][0]
            else:
                break
            return foundValue
        elif mode == "horizontal":
            if not clickValue['result'][1] is None:
                foundValue = clickValue['result'][1]
            else:
                break
            return foundValue


def id_func(specsfile):
    filename = os.path.basename(specsfile)
    print(filename)
    splitted_file = filename.split("_")
    print(splitted_file)
    starname = splitted_file[0]
    date = splitted_file[2]
    return starname, date


def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Run RV fitter from command line')
    parser.add_argument('--debug', action='store_true')
    # parser.add_argument('--threshold', action='store_true')
    # parser.add_argument('--asimov', action='store_true')
    parser.add_argument('--specfile_list', type=str, required=True)
    parser.add_argument('--line_list', type=str, required=True)
    # parser.add_argument('--seed', type=int, default=None, required=False)
    # parser.add_argument('--livetime', type=float, default=None)
    # parser.add_argument('--var1', type=str, required=True)
    # parser.add_argument('--val1', type=float, required=True)
    #  parser.add_argument('--var2', type=str, default=None, required=False)
    #  parser.add_argument('--val2', type=float, default=[None], required=False)
    #  parser.add_argument('--output', type=str, required=True)
    # parser.add_argument('--toy_size', type=int, default=500)
    return dict(vars(parser.parse_args()))


def continue_clipping(line):
    print('Do you want to continue clipping this line? (y/n)')
    global skip
    fig, ax = plt.subplots()
    line.plot_normed_spectrum(
        ax=ax,
        title_prefix=
        'Regions to be clipped\nDo you want to continue clipping this line? (y/n)\n'
    )
    line.plot_clips(ax=ax)
    plt.tight_layout()

    skip = {}
    skip['value'] = None
    fig.canvas.mpl_connect('key_press_event', lambda event: close_key(event))
    plt.show()
    return not skip['value']


def redo_entire_line(line):
    print('Do you want to re-do the normalizing and clipping? (y/n)')
    global skip
    fig, ax = plt.subplots()
    line.plot_clipped_spectrum(
        ax=ax,
        title_prefix=
        'Line profiles to be fitted\nDo you want to re-do the normalizing and clipping? (y/n)\n'
    )
    line.plot_clips(ax=ax)
    plt.tight_layout()

    skip = {}
    skip['value'] = None
    fig.canvas.mpl_connect('key_press_event', lambda event: close_key(event))
    plt.show()
    return not skip['value']


def main(args):
    parsed_args = parse_args(args)
    print(parsed_args['specfile_list'])
    line_list = parsed_args['line_list']
    with open(parsed_args['specfile_list'], 'r') as f:
        specsfilelist = f.read().splitlines()
    debug = parsed_args['debug']

    if debug:
        specsfilelist = specsfilelist[:1]

    myfitter = RVFitter.from_specsfilelist_name_flexi(
        specsfilelist_name=parsed_args['specfile_list'],
        id_func=id_func,
        line_list=line_list,
        debug=debug)

    for idx, rvobject in enumerate(myfitter.rvobjects):
        if idx != 0:
            rvobject.apply_selecting(standard_epoch=myfitter.rvobjects[0])
        for line in rvobject.lines:
            while True:
                # choosing lines
                if idx == 0:
                    skip = choose_line(rvObject=rvobject, line=line)
                    line.is_selected = not skip
                else:
                    skip = not line.is_selected
                if not skip:
                    line._clear()
                    # Normalizing
                    leftvalue_horizontal = get_clickValue(rvObject=rvobject,
                                                          line=line,
                                                          mode='horizontal',
                                                          position='left')
                    rightvalue_horizontal = get_clickValue(
                        rvObject=rvobject,
                        line=line,
                        mode='horizontal',
                        position='right',
                        previous_value=leftvalue_horizontal)

                    line.add_normed_spectrum(
                        angstrom=rvobject.angstrom,
                        flux=rvobject.flux,
                        error=rvobject.flux_errors,
                        leftValueNorm=leftvalue_horizontal,
                        rightValueNorm=rightvalue_horizontal)
                    # Clipping
                    while True:
                        leftvalue_vertical = get_clickValue(rvObject=rvobject,
                                                            line=line,
                                                            mode='vertical',
                                                            position='left')
                        rightvalue_vertical = get_clickValue(
                            rvObject=rvobject,
                            line=line,
                            mode='vertical',
                            position='right',
                            previous_value=leftvalue_vertical)

                        line.clip_spectrum(leftClip=leftvalue_vertical,
                                           rightClip=rightvalue_vertical)

                        keep_clipping = continue_clipping(line=line)
                        if not keep_clipping:
                            break

                else:
                    break

                repeat_process = redo_entire_line(line=line)
                if not repeat_process:
                    break

    myfitter.create_df()
    myfitter.save_df()
    print(myfitter.df)


if __name__ == "__main__":
    main(sys.argv[1:])
