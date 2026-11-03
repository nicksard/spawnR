package com.nicksard.appsleep;

import android.app.Activity;
import android.app.AlertDialog;
import android.app.TimePickerDialog;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.content.pm.ResolveInfo;
import android.graphics.drawable.Drawable;
import android.os.Bundle;
import android.text.Editable;
import android.text.InputType;
import android.text.TextWatcher;
import android.text.format.DateFormat;
import android.view.Menu;
import android.view.MenuItem;
import android.view.ViewGroup;
import android.widget.Button;
import android.widget.CheckBox;
import android.widget.EditText;
import android.widget.ImageView;
import android.widget.LinearLayout;
import android.widget.ScrollView;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.ToggleButton;

import java.text.Collator;
import java.time.DayOfWeek;
import java.time.format.TextStyle;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

public class EditScheduleActivity extends Activity {
    static final String EXTRA_ID = "id";
    private static final int MENU_SAVE = 1;
    private static final int MENU_DELETE = 2;

    private Schedule schedule;
    private boolean isNew;

    private EditText nameField;
    private Button startButton;
    private Button endButton;
    private TextView appsHeading;
    private LinearLayout appList;
    private final List<AppRow> appRows = new ArrayList<>();

    private static final class AppRow {
        String pkg;
        String label;
        Drawable icon;
        LinearLayout view;
    }

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        schedule = ScheduleStore.find(this, getIntent().getStringExtra(EXTRA_ID));
        isNew = schedule == null;
        if (isNew) schedule = new Schedule();
        setTitle(isNew ? "New schedule" : "Edit schedule");

        ScrollView scroll = new ScrollView(this);
        scroll.setFitsSystemWindows(true);
        LinearLayout root = Ui.column(this);
        int pad = Ui.dp(this, 16);
        root.setPadding(pad, pad, pad, Ui.dp(this, 48));
        scroll.addView(root);

        // Name
        root.addView(Ui.heading(this, "Name"));
        nameField = new EditText(this);
        nameField.setInputType(InputType.TYPE_CLASS_TEXT | InputType.TYPE_TEXT_FLAG_CAP_SENTENCES);
        nameField.setText(schedule.name);
        nameField.setSingleLine(true);
        root.addView(nameField);

        // Times
        root.addView(Ui.heading(this, "When"));
        LinearLayout times = Ui.row(this);
        startButton = new Button(this);
        startButton.setAllCaps(false);
        startButton.setOnClickListener(v -> pickTime(true));
        endButton = new Button(this);
        endButton.setAllCaps(false);
        endButton.setOnClickListener(v -> pickTime(false));
        times.addView(startButton, Ui.weight1());
        times.addView(endButton, Ui.weight1());
        root.addView(times);
        updateTimeButtons();
        root.addView(Ui.secondary(this,
                "If the end is earlier than the start, the window runs overnight. "
                + "Same start and end means all day."));

        // Days
        root.addView(Ui.heading(this, "Days"));
        LinearLayout days = Ui.row(this);
        for (DayOfWeek d : DayOfWeek.values()) {
            ToggleButton tb = new ToggleButton(this);
            String label = d.getDisplayName(TextStyle.NARROW, Locale.getDefault());
            tb.setText(label);
            tb.setTextOn(label);
            tb.setTextOff(label);
            tb.setChecked(schedule.isDayOn(d));
            tb.setMinWidth(0);
            tb.setMinimumWidth(0);
            tb.setOnCheckedChangeListener((b, on) -> schedule.setDayOn(d, on));
            days.addView(tb, Ui.weight1());
        }
        root.addView(days);

        // Notifications
        CheckBox silence = new CheckBox(this);
        silence.setText("Also clear their notifications while asleep");
        silence.setChecked(schedule.silenceNotifications);
        silence.setOnCheckedChangeListener((b, on) -> schedule.silenceNotifications = on);
        LinearLayout.LayoutParams slp = new LinearLayout.LayoutParams(
                ViewGroup.LayoutParams.MATCH_PARENT, ViewGroup.LayoutParams.WRAP_CONTENT);
        slp.topMargin = Ui.dp(this, 12);
        root.addView(silence, slp);

        // Apps
        appsHeading = Ui.heading(this, "Apps to put to sleep");
        root.addView(appsHeading);
        EditText search = new EditText(this);
        search.setHint("Search apps");
        search.setSingleLine(true);
        search.addTextChangedListener(new TextWatcher() {
            @Override public void beforeTextChanged(CharSequence s, int a, int b, int c) {}
            @Override public void onTextChanged(CharSequence s, int a, int b, int c) {}
            @Override public void afterTextChanged(Editable s) { filterApps(s.toString()); }
        });
        root.addView(search);
        appList = Ui.column(this);
        appList.addView(Ui.secondary(this, "Loading apps…"));
        root.addView(appList);

        Button save = new Button(this);
        save.setText("Save");
        save.setOnClickListener(v -> save());
        LinearLayout.LayoutParams blp = new LinearLayout.LayoutParams(
                ViewGroup.LayoutParams.MATCH_PARENT, ViewGroup.LayoutParams.WRAP_CONTENT);
        blp.topMargin = Ui.dp(this, 16);
        root.addView(save, blp);

        setContentView(scroll);
        loadApps();
    }

    private void updateTimeButtons() {
        startButton.setText("Sleeps at\n" + Schedule.formatMinute(schedule.startMinute));
        endButton.setText("Wakes at\n" + Schedule.formatMinute(schedule.endMinute));
    }

    private void updateAppsHeading() {
        int n = schedule.packages.size();
        appsHeading.setText("Apps to put to sleep" + (n > 0 ? " (" + n + " selected)" : ""));
    }

    private void pickTime(boolean start) {
        int m = start ? schedule.startMinute : schedule.endMinute;
        new TimePickerDialog(this, (view, hour, minute) -> {
            if (start) schedule.startMinute = hour * 60 + minute;
            else schedule.endMinute = hour * 60 + minute;
            updateTimeButtons();
        }, m / 60, m % 60, DateFormat.is24HourFormat(this)).show();
    }

    private void loadApps() {
        new Thread(() -> {
            PackageManager pm = getPackageManager();
            Intent launcher = new Intent(Intent.ACTION_MAIN).addCategory(Intent.CATEGORY_LAUNCHER);
            List<ResolveInfo> infos = pm.queryIntentActivities(launcher, 0);
            Set<String> seen = new HashSet<>();
            List<AppRow> rows = new ArrayList<>();
            for (ResolveInfo ri : infos) {
                String pkg = ri.activityInfo.packageName;
                if (pkg.equals(getPackageName()) || !seen.add(pkg)) continue;
                AppRow r = new AppRow();
                r.pkg = pkg;
                r.label = ri.loadLabel(pm).toString();
                r.icon = ri.loadIcon(pm);
                rows.add(r);
            }
            // Apps already chosen first, then alphabetical.
            Collator collator = Collator.getInstance();
            rows.sort((a, b) -> {
                boolean sa = schedule.packages.contains(a.pkg);
                boolean sb = schedule.packages.contains(b.pkg);
                if (sa != sb) return sa ? -1 : 1;
                return collator.compare(a.label, b.label);
            });
            runOnUiThread(() -> showApps(rows));
        }).start();
    }

    private void showApps(List<AppRow> rows) {
        if (isFinishing() || isDestroyed()) return;
        appList.removeAllViews();
        appRows.clear();
        int iconSize = Ui.dp(this, 36);
        for (AppRow r : rows) {
            LinearLayout line = Ui.row(this);
            line.setPadding(0, Ui.dp(this, 6), 0, Ui.dp(this, 6));
            ImageView iv = new ImageView(this);
            iv.setImageDrawable(r.icon);
            LinearLayout.LayoutParams ilp = new LinearLayout.LayoutParams(iconSize, iconSize);
            ilp.rightMargin = Ui.dp(this, 12);
            line.addView(iv, ilp);
            CheckBox cb = new CheckBox(this);
            cb.setText(r.label);
            cb.setChecked(schedule.packages.contains(r.pkg));
            cb.setOnCheckedChangeListener((b, on) -> {
                if (on) schedule.packages.add(r.pkg);
                else schedule.packages.remove(r.pkg);
                updateAppsHeading();
            });
            line.addView(cb, Ui.weight1());
            line.setOnClickListener(v -> cb.toggle());
            r.view = line;
            appRows.add(r);
            appList.addView(line);
        }
        updateAppsHeading();
    }

    private void filterApps(String query) {
        String q = query.trim().toLowerCase(Locale.getDefault());
        for (AppRow r : appRows) {
            boolean match = q.isEmpty()
                    || r.label.toLowerCase(Locale.getDefault()).contains(q)
                    || r.pkg.toLowerCase(Locale.ROOT).contains(q);
            Ui.gone(r.view, !match);
        }
    }

    private void save() {
        String name = nameField.getText().toString().trim();
        schedule.name = name.isEmpty() ? "Sleep schedule" : name;
        if (schedule.daysMask == 0) {
            Toast.makeText(this, "Pick at least one day.", Toast.LENGTH_SHORT).show();
            return;
        }
        if (schedule.packages.isEmpty()) {
            Toast.makeText(this, "Pick at least one app to put to sleep.", Toast.LENGTH_SHORT).show();
            return;
        }
        ScheduleStore.save(this, schedule);
        finish();
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        menu.add(0, MENU_SAVE, 0, "Save").setShowAsAction(MenuItem.SHOW_AS_ACTION_ALWAYS);
        if (!isNew) {
            menu.add(0, MENU_DELETE, 1, "Delete").setShowAsAction(MenuItem.SHOW_AS_ACTION_NEVER);
        }
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        if (item.getItemId() == MENU_SAVE) {
            save();
            return true;
        }
        if (item.getItemId() == MENU_DELETE) {
            new AlertDialog.Builder(this)
                    .setTitle("Delete \"" + schedule.name + "\"?")
                    .setPositiveButton("Delete", (d, w) -> {
                        ScheduleStore.delete(this, schedule.id);
                        finish();
                    })
                    .setNegativeButton("Cancel", null)
                    .show();
            return true;
        }
        return super.onOptionsItemSelected(item);
    }
}
